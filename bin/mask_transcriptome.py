#!/usr/bin/env python3
"""
mask_transcriptome.py
─────────────────────────────────────────────────────────────────────────────
HeteroPhy Pipeline — Step 5: Transcriptome Masking

Replaces nucleotides at positions with insufficient read coverage with 'N'.

RATIONALE:
    In phylogenetic likelihood calculations:
    - 'N' (unknown state) = treated as all states simultaneously → no bias
    - Wrong base          = actively misleads topology inference

    This is especially critical for rapidly radiating species where the
    signal resides in a small number of differences between species.

Author  : HeteroPhy Pipeline
Version : 1.0.0
─────────────────────────────────────────────────────────────────────────────
"""

import argparse
import json
import sys
from pathlib import Path
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(
        description="Mask low-coverage nucleotide positions in transcriptome FASTA",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument("--fasta",      required=True, help="Input transcriptome FASTA")
    parser.add_argument("--coverage",   required=True, help="Samtools depth TSV (contig, pos, depth)")
    parser.add_argument("--min_depth",  type=int, default=5, help="Minimum depth to keep a position")
    parser.add_argument("--min_length", type=int, default=200, help="Minimum contig length after masking")
    parser.add_argument("--sample",     required=True, help="Sample ID for logging")
    parser.add_argument("--output",     required=True, help="Output masked FASTA")
    parser.add_argument("--stats",      required=True, help="Output JSON statistics file")
    return parser.parse_args()


def read_fasta(path):
    """Parse FASTA into ordered dict {header: sequence}."""
    sequences = {}
    order = []
    current_header = None
    current_seq = []

    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if current_header is not None:
                    sequences[current_header] = "".join(current_seq)
                current_header = line[1:].split()[0]   # use sequence ID only
                current_seq = []
                order.append(current_header)
            else:
                current_seq.append(line.upper())

    if current_header is not None:
        sequences[current_header] = "".join(current_seq)

    return sequences, order


def read_coverage(path):
    """
    Parse samtools depth output into dict: {contig: {position: depth}}
    Samtools depth: col1=contig, col2=pos(1-based), col3=depth
    """
    coverage = defaultdict(dict)
    with open(path) as f:
        for line in f:
            parts = line.rstrip().split("\t")
            if len(parts) < 3:
                continue
            contig, pos, depth = parts[0], int(parts[1]), int(parts[2])
            coverage[contig][pos] = depth
    return coverage


def mask_sequence(seq, coverage_dict, min_depth):
    """
    Replace positions with depth < min_depth with 'N'.
    Returns: masked sequence, number of positions masked, fraction masked.
    """
    seq_list = list(seq)
    n_masked = 0

    for i, base in enumerate(seq_list):
        pos = i + 1  # 1-based
        depth = coverage_dict.get(pos, 0)
        if depth < min_depth:
            seq_list[i] = "N"
            n_masked += 1

    return "".join(seq_list), n_masked


def calculate_n_fraction(seq):
    """Return fraction of sequence that is N or ambiguous."""
    n_count = sum(1 for b in seq if b == "N")
    return n_count / len(seq) if seq else 0.0


def main():
    args = parse_args()

    print(f"[HeteroPhy] Masking transcriptome for sample: {args.sample}")
    print(f"[HeteroPhy] Minimum coverage depth: {args.min_depth}")
    print(f"[HeteroPhy] Minimum contig length post-masking: {args.min_length}")

    # ── Load inputs ───────────────────────────────────────────────────────────
    sequences, order = read_fasta(args.fasta)
    coverage = read_coverage(args.coverage)

    print(f"[HeteroPhy] Loaded {len(sequences):,} contigs from FASTA")
    print(f"[HeteroPhy] Loaded coverage data for {len(coverage):,} contigs")

    # ── Mask each sequence ────────────────────────────────────────────────────
    masked_sequences = {}
    stats = {
        "sample": args.sample,
        "total_contigs": len(sequences),
        "retained_contigs": 0,
        "filtered_contigs": 0,
        "total_bases": 0,
        "masked_bases": 0,
        "fraction_masked": 0.0,
        "contigs_no_coverage": 0,
        "min_depth_threshold": args.min_depth,
        "min_length_threshold": args.min_length,
    }

    for contig_id in order:
        seq = sequences[contig_id]
        cov_dict = coverage.get(contig_id, {})

        if not cov_dict:
            stats["contigs_no_coverage"] += 1
            # Entire contig has no coverage → all N
            masked_seq = "N" * len(seq)
            n_masked = len(seq)
        else:
            masked_seq, n_masked = mask_sequence(seq, cov_dict, args.min_depth)

        stats["total_bases"] += len(seq)
        stats["masked_bases"] += n_masked

        # Filter out contigs that are too short or completely masked
        non_n = len(masked_seq) - masked_seq.count("N")
        if non_n >= args.min_length:
            masked_sequences[contig_id] = masked_seq
            stats["retained_contigs"] += 1
        else:
            stats["filtered_contigs"] += 1

    # ── Summary stats ─────────────────────────────────────────────────────────
    if stats["total_bases"] > 0:
        stats["fraction_masked"] = round(stats["masked_bases"] / stats["total_bases"], 4)

    print(f"[HeteroPhy] Retained: {stats['retained_contigs']:,} contigs")
    print(f"[HeteroPhy] Filtered: {stats['filtered_contigs']:,} contigs (below min_length)")
    print(f"[HeteroPhy] Bases masked: {stats['masked_bases']:,} / {stats['total_bases']:,} "
          f"({stats['fraction_masked']*100:.1f}%)")

    # ── Write masked FASTA ────────────────────────────────────────────────────
    with open(args.output, "w") as out:
        for contig_id, seq in masked_sequences.items():
            out.write(f">{contig_id}\n")
            # Wrap at 80 characters
            for i in range(0, len(seq), 80):
                out.write(seq[i:i+80] + "\n")

    # ── Write stats JSON ──────────────────────────────────────────────────────
    with open(args.stats, "w") as out:
        json.dump(stats, out, indent=2)

    print(f"[HeteroPhy] Masked FASTA written to: {args.output}")
    print(f"[HeteroPhy] Stats written to: {args.stats}")


if __name__ == "__main__":
    main()

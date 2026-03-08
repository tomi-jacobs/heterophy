#!/usr/bin/env python3
"""
apply_heterozygosity_strategy.py
─────────────────────────────────────────────────────────────────────────────
HeteroPhy Pipeline — Steps 8a/8b: Heterozygosity Strategies

Applies one of two strategies to handle heterozygous SNP sites:

  STRATEGY 1 — COLLAPSE:
      At heterozygous positions, the reference allele (or major allele by
      read depth) is used. Simplifies data for traditional phylogenetic
      models. Information loss, but avoids introducing non-existent
      genotypes into the inference.

  STRATEGY 2 — RETAIN:
      Heterozygous positions are encoded using IUPAC ambiguity codes.
      Preserves allelic variation that may reflect recent divergence.
      Risk: may introduce phantom genotypes not present in either true
      haplotype.

IUPAC ambiguity code table used:
    R = A/G    Y = C/T    S = G/C
    W = A/T    K = G/T    M = A/C
    B = C/G/T  D = A/G/T  H = A/C/T
    V = A/C/G  N = any

Author  : HeteroPhy Pipeline
Version : 1.0.0
─────────────────────────────────────────────────────────────────────────────
"""

import argparse
import gzip
import json
import sys
from pathlib import Path
from collections import defaultdict


# ── IUPAC ambiguity code lookup ───────────────────────────────────────────────
IUPAC = {
    frozenset(["A", "G"]): "R",
    frozenset(["C", "T"]): "Y",
    frozenset(["G", "C"]): "S",
    frozenset(["A", "T"]): "W",
    frozenset(["G", "T"]): "K",
    frozenset(["A", "C"]): "M",
    frozenset(["C", "G", "T"]): "B",
    frozenset(["A", "G", "T"]): "D",
    frozenset(["A", "C", "T"]): "H",
    frozenset(["A", "C", "G"]): "V",
    frozenset(["A"]): "A",
    frozenset(["C"]): "C",
    frozenset(["G"]): "G",
    frozenset(["T"]): "T",
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Apply heterozygosity strategy (collapse or retain) to masked FASTA",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument("--fasta",    required=True, help="Masked transcriptome FASTA")
    parser.add_argument("--vcf",      required=True, help="Filtered VCF (gzipped)")
    parser.add_argument("--strategy", required=True, choices=["collapse", "retain"],
                        help="Heterozygosity handling strategy")
    parser.add_argument("--sample",   required=True, help="Sample ID")
    parser.add_argument("--output",   required=True, help="Output FASTA")
    return parser.parse_args()


def read_fasta(path):
    """Parse FASTA into dict {header: sequence}."""
    sequences = {}
    order = []
    current_header = None
    current_seq = []

    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if current_header is not None:
                    sequences[current_header] = list("".join(current_seq))
                current_header = line[1:].split()[0]
                current_seq = []
                order.append(current_header)
            else:
                current_seq.append(line.upper())

    if current_header is not None:
        sequences[current_header] = list("".join(current_seq))

    return sequences, order


def parse_vcf_heterozygous(vcf_path, sample_col=9):
    """
    Parse VCF for heterozygous sites.
    Returns: dict {contig: {pos(0-based): (ref, alt, gt_type)}}
    gt_type: 'het' or 'hom_alt'
    """
    het_sites = defaultdict(dict)
    opener = gzip.open if vcf_path.endswith(".gz") else open

    with opener(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip().split("\t")
            if len(parts) < 10:
                continue

            chrom  = parts[0]
            pos    = int(parts[1]) - 1   # Convert to 0-based
            ref    = parts[3].upper()
            alts   = parts[4].upper().split(",")
            fmt    = parts[8].split(":")
            sample = parts[sample_col].split(":")

            # Parse genotype
            if "GT" not in fmt:
                continue
            gt_idx = fmt.index("GT")
            if gt_idx >= len(sample):
                continue

            gt = sample[gt_idx].replace("|", "/")
            if "." in gt:
                continue   # Missing genotype

            alleles = gt.split("/")
            if len(alleles) != 2:
                continue

            # Get allele sequences
            all_alleles = [ref] + alts
            try:
                a1 = all_alleles[int(alleles[0])]
                a2 = all_alleles[int(alleles[1])]
            except (IndexError, ValueError):
                continue

            # Only process SNPs (single base)
            if len(ref) != 1 or len(a1) != 1 or len(a2) != 1:
                continue

            if alleles[0] != alleles[1]:
                het_sites[chrom][pos] = (ref, a1, a2, "het")
            elif alleles[0] != "0":
                het_sites[chrom][pos] = (ref, a1, a2, "hom_alt")

    return het_sites


def apply_collapse(sequences, order, het_sites):
    """
    COLLAPSE STRATEGY:
    At heterozygous sites, use the reference allele.
    At homozygous alt sites, use the alt allele.
    """
    stats = {"het_sites_collapsed": 0, "hom_alt_applied": 0}

    for contig in order:
        if contig not in sequences:
            continue
        sites = het_sites.get(contig, {})
        for pos, (ref, a1, a2, gt_type) in sites.items():
            if pos >= len(sequences[contig]):
                continue
            if gt_type == "het":
                sequences[contig][pos] = ref   # Use reference at het sites
                stats["het_sites_collapsed"] += 1
            elif gt_type == "hom_alt":
                sequences[contig][pos] = a1    # Apply homozygous alt
                stats["hom_alt_applied"] += 1

    return sequences, stats


def apply_retain(sequences, order, het_sites):
    """
    RETAIN STRATEGY:
    At heterozygous sites, encode both alleles using IUPAC ambiguity codes.
    At homozygous alt sites, apply the alt allele.
    """
    stats = {"het_sites_encoded": 0, "iupac_unknown": 0, "hom_alt_applied": 0}

    for contig in order:
        if contig not in sequences:
            continue
        sites = het_sites.get(contig, {})
        for pos, (ref, a1, a2, gt_type) in sites.items():
            if pos >= len(sequences[contig]):
                continue
            if gt_type == "het":
                key = frozenset([a1, a2])
                iupac_code = IUPAC.get(key, "N")   # Fall back to N if not found
                if iupac_code == "N":
                    stats["iupac_unknown"] += 1
                sequences[contig][pos] = iupac_code
                stats["het_sites_encoded"] += 1
            elif gt_type == "hom_alt":
                sequences[contig][pos] = a1
                stats["hom_alt_applied"] += 1

    return sequences, stats


def write_fasta(sequences, order, output_path):
    with open(output_path, "w") as out:
        for contig_id in order:
            if contig_id not in sequences:
                continue
            seq = "".join(sequences[contig_id])
            out.write(f">{contig_id}\n")
            for i in range(0, len(seq), 80):
                out.write(seq[i:i+80] + "\n")


def main():
    args = parse_args()

    print(f"[HeteroPhy] Applying strategy: {args.strategy.upper()} for sample: {args.sample}")

    # ── Load inputs ───────────────────────────────────────────────────────────
    sequences, order = read_fasta(args.fasta)
    het_sites = parse_vcf_heterozygous(args.vcf)

    total_het = sum(len(v) for v in het_sites.values())
    print(f"[HeteroPhy] Loaded {len(sequences):,} contigs")
    print(f"[HeteroPhy] Found {total_het:,} variant sites across "
          f"{len(het_sites):,} contigs")

    # ── Apply strategy ────────────────────────────────────────────────────────
    if args.strategy == "collapse":
        sequences, stats = apply_collapse(sequences, order, het_sites)
        print(f"[HeteroPhy] Collapsed {stats['het_sites_collapsed']:,} het sites "
              f"(retained reference allele)")
        print(f"[HeteroPhy] Applied {stats['hom_alt_applied']:,} homozygous alt alleles")
    else:
        sequences, stats = apply_retain(sequences, order, het_sites)
        print(f"[HeteroPhy] Encoded {stats['het_sites_encoded']:,} het sites "
              f"with IUPAC ambiguity codes")
        if stats["iupac_unknown"] > 0:
            print(f"[HeteroPhy] WARNING: {stats['iupac_unknown']:,} sites fell back to N "
                  f"(multi-allelic or unrecognized)")
        print(f"[HeteroPhy] Applied {stats['hom_alt_applied']:,} homozygous alt alleles")

    # ── Write output ──────────────────────────────────────────────────────────
    write_fasta(sequences, order, args.output)
    print(f"[HeteroPhy] Output written to: {args.output}")


if __name__ == "__main__":
    main()

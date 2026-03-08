#!/usr/bin/env python3
"""
coverage_assessment.py
─────────────────────
HeteroPhy — Step 4: Nucleotide-level coverage assessment.
Parses samtools depth output and produces a BED file of low-coverage regions.
"""

import argparse
import json
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(description="Assess per-nucleotide coverage")
    parser.add_argument("--coverage",  required=True, help="Samtools depth TSV")
    parser.add_argument("--min_depth", type=int, default=5)
    parser.add_argument("--sample",    required=True)
    parser.add_argument("--bed_out",   required=True)
    parser.add_argument("--stats_out", required=True)
    return parser.parse_args()


def main():
    args = parse_args()

    # Track per-contig coverage
    contig_positions = defaultdict(list)
    total_positions = 0
    low_cov_positions = 0

    with open(args.coverage) as f:
        for line in f:
            parts = line.rstrip().split("\t")
            if len(parts) < 3:
                continue
            contig, pos, depth = parts[0], int(parts[1]), int(parts[2])
            contig_positions[contig].append((pos, depth))
            total_positions += 1
            if depth < args.min_depth:
                low_cov_positions += 1

    # Write BED file of low-coverage regions (merge consecutive positions)
    with open(args.bed_out, "w") as bed:
        for contig, positions in contig_positions.items():
            positions.sort(key=lambda x: x[0])
            region_start = None
            region_end = None

            for pos, depth in positions:
                if depth < args.min_depth:
                    if region_start is None:
                        region_start = pos - 1   # BED is 0-based
                    region_end = pos
                else:
                    if region_start is not None:
                        bed.write(f"{contig}\t{region_start}\t{region_end}\n")
                        region_start = None
                        region_end = None

            if region_start is not None:
                bed.write(f"{contig}\t{region_start}\t{region_end}\n")

    # Compute per-contig stats
    contig_stats = {}
    for contig, positions in contig_positions.items():
        depths = [d for _, d in positions]
        if depths:
            contig_stats[contig] = {
                "length": len(depths),
                "mean_depth": round(sum(depths) / len(depths), 2),
                "min_depth": min(depths),
                "max_depth": max(depths),
                "low_cov_bases": sum(1 for d in depths if d < args.min_depth),
            }

    stats = {
        "sample": args.sample,
        "total_positions": total_positions,
        "low_coverage_positions": low_cov_positions,
        "fraction_low_coverage": round(low_cov_positions / total_positions, 4)
            if total_positions > 0 else 0.0,
        "n_contigs": len(contig_positions),
        "min_depth_threshold": args.min_depth,
        "per_contig": contig_stats,
    }

    with open(args.stats_out, "w") as f:
        json.dump(stats, f, indent=2)

    print(f"[HeteroPhy] Coverage assessment complete for {args.sample}")
    print(f"[HeteroPhy] {low_cov_positions:,}/{total_positions:,} positions "
          f"below min depth ({args.min_depth}x)")


if __name__ == "__main__":
    main()

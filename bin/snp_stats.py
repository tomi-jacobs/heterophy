#!/usr/bin/env python3
"""
snp_stats.py
────────────
HeteroPhy — Summarise per-sample SNP statistics from a filtered VCF.
"""

import argparse
import gzip
import json
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(description="Summarise SNP statistics from VCF")
    parser.add_argument("--vcf",    required=True, help="Filtered VCF (gzipped)")
    parser.add_argument("--sample", required=True, help="Sample ID")
    parser.add_argument("--output", required=True, help="Output JSON stats file")
    return parser.parse_args()


def main():
    args = parse_args()

    stats = {
        "sample": args.sample,
        "total_variants": 0,
        "snps": 0,
        "indels": 0,
        "heterozygous": 0,
        "homozygous_alt": 0,
        "homozygous_ref": 0,
        "transition_count": 0,
        "transversion_count": 0,
        "ts_tv_ratio": None,
        "per_contig_snp_density": {},
    }

    transitions = {frozenset(["A", "G"]), frozenset(["C", "T"])}
    contig_snp_counts = defaultdict(int)
    contig_lengths = defaultdict(int)

    opener = gzip.open if args.vcf.endswith(".gz") else open

    with opener(args.vcf, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue

            parts = line.rstrip().split("\t")
            if len(parts) < 10:
                continue

            chrom = parts[0]
            ref   = parts[3].upper()
            alts  = parts[4].upper().split(",")
            fmt   = parts[8].split(":")
            samp  = parts[9].split(":")

            stats["total_variants"] += 1

            # SNP vs indel
            is_snp = len(ref) == 1 and all(len(a) == 1 for a in alts if a != "*")
            if is_snp:
                stats["snps"] += 1
                contig_snp_counts[chrom] += 1

                # Transition / transversion
                for alt in alts:
                    if len(alt) == 1:
                        pair = frozenset([ref, alt])
                        if pair in transitions:
                            stats["transition_count"] += 1
                        else:
                            stats["transversion_count"] += 1
            else:
                stats["indels"] += 1

            # Genotype classification
            if "GT" in fmt:
                gt_idx = fmt.index("GT")
                if gt_idx < len(samp):
                    gt = samp[gt_idx].replace("|", "/")
                    alleles = gt.split("/")
                    if len(alleles) == 2 and "." not in gt:
                        if alleles[0] == alleles[1] == "0":
                            stats["homozygous_ref"] += 1
                        elif alleles[0] == alleles[1]:
                            stats["homozygous_alt"] += 1
                        else:
                            stats["heterozygous"] += 1

    # Ts/Tv ratio
    if stats["transversion_count"] > 0:
        stats["ts_tv_ratio"] = round(
            stats["transition_count"] / stats["transversion_count"], 3
        )

    # Per-contig density (if we had contig lengths — simplified here)
    stats["per_contig_snp_counts"] = dict(contig_snp_counts)
    stats["n_contigs_with_snps"] = len(contig_snp_counts)

    with open(args.output, "w") as f:
        json.dump(stats, f, indent=2)

    print(f"[HeteroPhy] SNP stats for {args.sample}: "
          f"{stats['snps']:,} SNPs, {stats['heterozygous']:,} het, "
          f"Ts/Tv={stats['ts_tv_ratio']}")


if __name__ == "__main__":
    main()

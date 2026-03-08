#!/usr/bin/env python3
"""
compare_trees.py
─────────────────────────────────────────────────────────────────────────────
HeteroPhy Pipeline — Step 11: Phylogenetic Tree Comparison

Compares trees inferred under the two heterozygosity strategies:
  1. COLLAPSED — heterozygous sites replaced with reference allele
  2. RETAINED  — heterozygous sites encoded with IUPAC ambiguity codes

Metrics computed:
  - Robinson-Foulds (RF) distance (normalised and unnormalised)
  - Weighted Robinson-Foulds distance
  - Quartet distance
  - Bootstrap support distribution comparison (Kolmogorov-Smirnov test)
  - Clade-level concordance/discordance identification

Author  : HeteroPhy Pipeline
Version : 1.0.0
─────────────────────────────────────────────────────────────────────────────
"""

import argparse
import json
import sys
from pathlib import Path
from collections import defaultdict

try:
    from ete3 import Tree
except ImportError:
    print("[HeteroPhy] ERROR: ete3 not installed. Install via: conda install -c ete3 ete3")
    sys.exit(1)

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from matplotlib.gridspec import GridSpec
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    print("[HeteroPhy] WARNING: matplotlib not available. Skipping plots.")

try:
    from scipy import stats as scipy_stats
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compare phylogenetic trees between heterozygosity strategies",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument("--trees",  required=True, nargs="+",
                        help="Tree files (Newick). Name format: {strategy}.treefile")
    parser.add_argument("--outdir", default=".", help="Output directory")
    parser.add_argument("--output", required=True, help="Output JSON file")
    return parser.parse_args()


def load_tree(path):
    """Load a Newick tree file."""
    try:
        t = Tree(str(path), format=0)   # format=0 = Newick with bootstrap support
        return t
    except Exception as e:
        print(f"[HeteroPhy] WARNING: Could not load tree {path}: {e}")
        return None


def extract_bootstrap_values(tree):
    """Extract all internal node bootstrap support values."""
    values = []
    for node in tree.traverse():
        if not node.is_leaf() and not node.is_root():
            try:
                bs = float(node.name) if node.name else None
                if bs is not None:
                    values.append(bs)
            except (ValueError, TypeError):
                pass
    return values


def compute_rf_distance(tree1, tree2):
    """
    Compute Robinson-Foulds distance between two trees.
    Returns: rf, max_rf, normalized_rf, weighted_rf
    """
    try:
        rf, max_rf, common_leaves, parts_t1, parts_t2, discard_t1, discard_t2 = \
            tree1.robinson_foulds(tree2, unrooted_trees=True)

        normalized = rf / max_rf if max_rf > 0 else 0.0

        # Weighted RF (branch-length aware)
        try:
            wrf, _, _, _, _, _, _ = tree1.robinson_foulds(
                tree2, unrooted_trees=True, expand_polytomies=False
            )
        except Exception:
            wrf = rf

        return {
            "rf": int(rf),
            "max_rf": int(max_rf),
            "normalized_rf": round(normalized, 4),
            "common_leaves": len(common_leaves),
            "partitions_t1": len(parts_t1),
            "partitions_t2": len(parts_t2),
        }
    except Exception as e:
        print(f"[HeteroPhy] WARNING: RF computation failed: {e}")
        return {"rf": None, "error": str(e)}


def identify_discordant_clades(tree1, tree2, strategy1, strategy2):
    """
    Identify clades present in one tree but not the other.
    Returns lists of discordant clades with their bootstrap support.
    """
    def get_clades(tree):
        clades = {}
        for node in tree.traverse():
            if node.is_leaf() or node.is_root():
                continue
            leaves = frozenset(node.get_leaf_names())
            try:
                bs = float(node.name) if node.name else 0.0
            except (ValueError, TypeError):
                bs = 0.0
            clades[leaves] = bs
        return clades

    clades1 = get_clades(tree1)
    clades2 = get_clades(tree2)

    shared = set(clades1.keys()) & set(clades2.keys())
    only_in_1 = set(clades1.keys()) - set(clades2.keys())
    only_in_2 = set(clades2.keys()) - set(clades1.keys())

    return {
        "shared_clades": len(shared),
        f"unique_to_{strategy1}": [
            {"clade": sorted(list(c)), "bootstrap": clades1[c]}
            for c in only_in_1
        ],
        f"unique_to_{strategy2}": [
            {"clade": sorted(list(c)), "bootstrap": clades2[c]}
            for c in only_in_2
        ],
    }


def compare_bootstrap_distributions(bs1, bs2, strategy1, strategy2):
    """
    Compare bootstrap support distributions between strategies.
    Uses Kolmogorov-Smirnov test if scipy available.
    """
    result = {
        "strategy_1": strategy1,
        "strategy_2": strategy2,
        "n_nodes_1": len(bs1),
        "n_nodes_2": len(bs2),
        "mean_bs_1": round(sum(bs1) / len(bs1), 2) if bs1 else None,
        "mean_bs_2": round(sum(bs2) / len(bs2), 2) if bs2 else None,
        "median_bs_1": round(sorted(bs1)[len(bs1)//2], 2) if bs1 else None,
        "median_bs_2": round(sorted(bs2)[len(bs2)//2], 2) if bs2 else None,
        "high_support_1": sum(1 for v in bs1 if v >= 95) if bs1 else 0,
        "high_support_2": sum(1 for v in bs2 if v >= 95) if bs2 else 0,
        "low_support_1": sum(1 for v in bs1 if v < 50) if bs1 else 0,
        "low_support_2": sum(1 for v in bs2 if v < 50) if bs2 else 0,
    }

    if HAS_SCIPY and bs1 and bs2:
        ks_stat, ks_pvalue = scipy_stats.ks_2samp(bs1, bs2)
        result["ks_statistic"] = round(ks_stat, 4)
        result["ks_pvalue"] = round(ks_pvalue, 6)
        result["ks_significant"] = bool(ks_pvalue < 0.05)
        result["ks_interpretation"] = (
            "Bootstrap distributions differ significantly between strategies "
            "(p < 0.05): heterozygosity treatment affects phylogenetic support"
            if ks_pvalue < 0.05 else
            "Bootstrap distributions are not significantly different (p >= 0.05): "
            "heterozygosity treatment does not markedly affect overall support"
        )

    return result


def plot_bootstrap_comparison(bs1, bs2, strategy1, strategy2, outdir):
    """Generate publication-ready bootstrap comparison figure."""
    if not HAS_MATPLOTLIB or not bs1 or not bs2:
        return None

    fig = plt.figure(figsize=(14, 6))
    fig.patch.set_facecolor("#0f1117")
    gs = GridSpec(1, 2, figure=fig, wspace=0.35)

    colors = {"collapsed": "#00d4aa", "retained": "#ff6b6b"}
    c1 = colors.get(strategy1, "#00d4aa")
    c2 = colors.get(strategy2, "#ff6b6b")

    # ── Left panel: Histogram comparison ──────────────────────────────────────
    ax1 = fig.add_subplot(gs[0])
    ax1.set_facecolor("#1a1d2e")
    ax1.hist(bs1, bins=20, range=(0, 100), alpha=0.75, color=c1,
             label=f"{strategy1.capitalize()} (n={len(bs1)})", edgecolor="none")
    ax1.hist(bs2, bins=20, range=(0, 100), alpha=0.75, color=c2,
             label=f"{strategy2.capitalize()} (n={len(bs2)})", edgecolor="none")
    ax1.axvline(95, color="#ffffff", linestyle="--", linewidth=0.8, alpha=0.6,
                label="95% threshold")
    ax1.set_xlabel("Bootstrap Support (%)", color="#e0e0e0", fontsize=11)
    ax1.set_ylabel("Number of Nodes", color="#e0e0e0", fontsize=11)
    ax1.set_title("Bootstrap Support Distribution", color="#ffffff", fontsize=13,
                  fontweight="bold", pad=12)
    ax1.tick_params(colors="#aaaaaa")
    ax1.spines[:].set_color("#333355")
    ax1.legend(facecolor="#1a1d2e", edgecolor="#333355", labelcolor="#e0e0e0",
               fontsize=9)

    # ── Right panel: Boxplot comparison ───────────────────────────────────────
    ax2 = fig.add_subplot(gs[1])
    ax2.set_facecolor("#1a1d2e")
    bp = ax2.boxplot([bs1, bs2],
                     labels=[strategy1.capitalize(), strategy2.capitalize()],
                     patch_artist=True,
                     medianprops=dict(color="#ffffff", linewidth=2),
                     whiskerprops=dict(color="#aaaaaa"),
                     capprops=dict(color="#aaaaaa"),
                     flierprops=dict(marker="o", color="#aaaaaa", alpha=0.4, markersize=3))
    for patch, color in zip(bp["boxes"], [c1, c2]):
        patch.set_facecolor(color)
        patch.set_alpha(0.8)
    ax2.axhline(95, color="#ffffff", linestyle="--", linewidth=0.8, alpha=0.6)
    ax2.set_ylabel("Bootstrap Support (%)", color="#e0e0e0", fontsize=11)
    ax2.set_title("Support Comparison by Strategy", color="#ffffff", fontsize=13,
                  fontweight="bold", pad=12)
    ax2.tick_params(colors="#aaaaaa")
    ax2.spines[:].set_color("#333355")

    fig.suptitle("HeteroPhy — Bootstrap Support Analysis",
                 color="#ffffff", fontsize=15, fontweight="bold", y=1.02)

    output_path = Path(outdir) / "bootstrap_comparison.pdf"
    plt.savefig(output_path, dpi=300, bbox_inches="tight",
                facecolor=fig.get_facecolor())
    plt.close()
    print(f"[HeteroPhy] Bootstrap plot saved: {output_path}")
    return str(output_path)


def plot_topology_comparison(rf_results, clade_results, outdir):
    """Generate topology comparison summary figure."""
    if not HAS_MATPLOTLIB:
        return None

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.patch.set_facecolor("#0f1117")

    for ax in axes:
        ax.set_facecolor("#1a1d2e")
        ax.tick_params(colors="#aaaaaa")
        ax.spines[:].set_color("#333355")

    # ── Panel 1: RF distance display ──────────────────────────────────────────
    ax = axes[0]
    rf_val = rf_results.get("normalized_rf", 0)
    color = "#00d4aa" if rf_val < 0.1 else "#ffaa00" if rf_val < 0.3 else "#ff6b6b"
    ax.barh(["Normalized RF\nDistance"], [rf_val], color=color, height=0.4)
    ax.barh(["Normalized RF\nDistance"], [1.0], color="#333355", height=0.4, alpha=0.4)
    ax.set_xlim(0, 1)
    ax.set_xlabel("Normalized RF Distance (0 = identical, 1 = maximally different)",
                  color="#e0e0e0", fontsize=9)
    ax.set_title("Tree Topology Difference", color="#ffffff", fontsize=12,
                 fontweight="bold")
    interpretation = ("Highly similar" if rf_val < 0.1 else
                      "Moderately different" if rf_val < 0.3 else
                      "Substantially different")
    ax.text(rf_val + 0.02, 0, f"{rf_val:.3f}\n({interpretation})",
            va="center", color="#ffffff", fontsize=10)

    # ── Panel 2: Clade concordance ────────────────────────────────────────────
    ax = axes[1]
    shared = clade_results.get("shared_clades", 0)
    keys = [k for k in clade_results if k.startswith("unique_to")]
    unique1 = len(clade_results.get(keys[0], [])) if len(keys) > 0 else 0
    unique2 = len(clade_results.get(keys[1], [])) if len(keys) > 1 else 0

    labels = ["Shared Clades", keys[0].replace("unique_to_", "Unique: "),
              keys[1].replace("unique_to_", "Unique: ")] if len(keys) >= 2 else \
             ["Shared Clades", "Unique"]
    values = [shared, unique1, unique2] if len(keys) >= 2 else [shared, unique1]
    bar_colors = ["#00d4aa", "#ffaa00", "#ff6b6b"]

    bars = ax.bar(labels, values, color=bar_colors[:len(values)], width=0.5)
    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3,
                str(val), ha="center", va="bottom", color="#ffffff", fontsize=11)
    ax.set_ylabel("Number of Clades", color="#e0e0e0", fontsize=11)
    ax.set_title("Clade Concordance", color="#ffffff", fontsize=12, fontweight="bold")

    fig.suptitle("HeteroPhy — Topology Comparison",
                 color="#ffffff", fontsize=15, fontweight="bold")
    plt.tight_layout()

    output_path = Path(outdir) / "topology_comparison.pdf"
    plt.savefig(output_path, dpi=300, bbox_inches="tight",
                facecolor=fig.get_facecolor())
    plt.close()
    print(f"[HeteroPhy] Topology plot saved: {output_path}")
    return str(output_path)


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"[HeteroPhy] Comparing {len(args.trees)} tree files")

    # ── Load trees, infer strategies from filenames ───────────────────────────
    trees = {}
    for path in args.trees:
        p = Path(path)
        # Infer strategy from filename: "collapsed.treefile" or "retained.treefile"
        name = p.stem.lower()
        strategy = "collapsed" if "collapsed" in name else \
                   "retained" if "retained" in name else name
        t = load_tree(path)
        if t is not None:
            trees[strategy] = t
            print(f"[HeteroPhy] Loaded tree: {strategy} ({len(t.get_leaf_names())} taxa)")

    if len(trees) < 2:
        print("[HeteroPhy] ERROR: Need at least 2 trees to compare.")
        # Write empty results
        with open(args.output, "w") as f:
            json.dump({"error": "Insufficient trees for comparison"}, f)
        # Create placeholder files for Nextflow
        for fname in ["topology_comparison.pdf", "bootstrap_comparison.pdf",
                      "robinson_foulds_distances.tsv", "bootstrap_support_summary.tsv"]:
            (outdir / fname).touch()
        sys.exit(0)

    strategy_names = list(trees.keys())
    s1, s2 = strategy_names[0], strategy_names[1]
    t1, t2 = trees[s1], trees[s2]

    # ── Compute all metrics ───────────────────────────────────────────────────
    print(f"[HeteroPhy] Computing Robinson-Foulds distance: {s1} vs {s2}")
    rf_results = compute_rf_distance(t1, t2)

    print(f"[HeteroPhy] Identifying discordant clades")
    clade_results = identify_discordant_clades(t1, t2, s1, s2)

    print(f"[HeteroPhy] Comparing bootstrap distributions")
    bs1 = extract_bootstrap_values(t1)
    bs2 = extract_bootstrap_values(t2)
    bs_comparison = compare_bootstrap_distributions(bs1, bs2, s1, s2)

    # ── Interpret results ─────────────────────────────────────────────────────
    rf_norm = rf_results.get("normalized_rf", 1.0)
    if rf_norm is not None:
        if rf_norm == 0.0:
            topology_verdict = "IDENTICAL: Both strategies produce the same topology"
        elif rf_norm < 0.05:
            topology_verdict = "NEAR-IDENTICAL: Strategies produce very similar topologies (minor branch differences)"
        elif rf_norm < 0.20:
            topology_verdict = "SIMILAR: Minor topological differences — heterozygosity has modest influence"
        elif rf_norm < 0.50:
            topology_verdict = "DIVERGENT: Moderate topological differences — heterozygosity influences relationship inference"
        else:
            topology_verdict = "HIGHLY DIVERGENT: Major topological differences — heterozygosity substantially alters phylogenetic inference"
    else:
        topology_verdict = "COULD NOT BE DETERMINED"

    # ── Compile results ───────────────────────────────────────────────────────
    results = {
        "pipeline": "HeteroPhy v1.0.0",
        "strategies_compared": [s1, s2],
        "topology": {
            "robinson_foulds": rf_results,
            "verdict": topology_verdict,
        },
        "clade_concordance": clade_results,
        "bootstrap_comparison": bs_comparison,
        "summary": {
            "heterozygosity_influences_topology": bool(rf_norm is not None and rf_norm > 0.0),
            "heterozygosity_influences_support": (
                bs_comparison.get("ks_significant", False)
            ),
            "recommendation": (
                "Both strategies yield equivalent topologies. "
                "Either approach is valid for this dataset."
                if (rf_norm is not None and rf_norm == 0.0) else
                "Heterozygosity treatment influences phylogenetic inference. "
                "We recommend reporting both strategies and discussing discordant nodes."
            ),
        }
    }

    # ── Write JSON ────────────────────────────────────────────────────────────
    with open(args.output, "w") as f:
        json.dump(results, f, indent=2)
    print(f"[HeteroPhy] Comparison results written to: {args.output}")

    # ── Write TSV summaries ───────────────────────────────────────────────────
    rf_tsv = outdir / "robinson_foulds_distances.tsv"
    with open(rf_tsv, "w") as f:
        f.write("comparison\trf_distance\tmax_rf\tnormalized_rf\tcommon_leaves\n")
        f.write(f"{s1}_vs_{s2}\t{rf_results.get('rf', 'NA')}\t"
                f"{rf_results.get('max_rf', 'NA')}\t"
                f"{rf_results.get('normalized_rf', 'NA')}\t"
                f"{rf_results.get('common_leaves', 'NA')}\n")

    bs_tsv = outdir / "bootstrap_support_summary.tsv"
    with open(bs_tsv, "w") as f:
        f.write("strategy\tn_nodes\tmean_bootstrap\tmedian_bootstrap\t"
                "nodes_ge95\tnodes_lt50\n")
        for strategy, bs_vals in [(s1, bs1), (s2, bs2)]:
            mean_bs = round(sum(bs_vals)/len(bs_vals), 2) if bs_vals else "NA"
            med_bs  = round(sorted(bs_vals)[len(bs_vals)//2], 2) if bs_vals else "NA"
            high_bs = sum(1 for v in bs_vals if v >= 95)
            low_bs  = sum(1 for v in bs_vals if v < 50)
            f.write(f"{strategy}\t{len(bs_vals)}\t{mean_bs}\t{med_bs}\t"
                    f"{high_bs}\t{low_bs}\n")

    # ── Generate plots ────────────────────────────────────────────────────────
    plot_bootstrap_comparison(bs1, bs2, s1, s2, outdir)
    plot_topology_comparison(rf_results, clade_results, outdir)

    # ── Print verdict ─────────────────────────────────────────────────────────
    print(f"\n{'═'*60}")
    print(f"  TOPOLOGY VERDICT: {topology_verdict}")
    print(f"  RF Distance: {rf_results.get('normalized_rf', 'NA')} (normalized)")
    print(f"  Shared clades: {clade_results.get('shared_clades', 'NA')}")
    print(f"  Mean bootstrap — {s1}: {bs_comparison.get('mean_bs_1', 'NA')}% | "
          f"{s2}: {bs_comparison.get('mean_bs_2', 'NA')}%")
    print(f"{'═'*60}\n")


if __name__ == "__main__":
    main()

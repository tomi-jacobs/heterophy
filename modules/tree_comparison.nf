process TREE_COMPARISON {
    tag "Comparing trees"
    label 'process_low'

    conda "envs/python.yml"
    container "heterophy/python:1.0.0"

    publishDir "${params.outdir}/comparison", mode: 'copy'

    input:
    path trees   // All .treefile outputs collected

    output:
    path "comparison_results.json",           emit: comparison_json
    path "topology_comparison.pdf",           emit: topology_plot
    path "bootstrap_comparison.pdf",          emit: bootstrap_plot
    path "robinson_foulds_distances.tsv",     emit: rf_distances
    path "bootstrap_support_summary.tsv",     emit: support_summary

    script:
    """
    compare_trees.py \\
        --trees      ${trees} \\
        --outdir     . \\
        --output     comparison_results.json
    """
}

/*
========================================================================================
    MODULE: generate_report — HTML report summarising the full analysis
========================================================================================
*/

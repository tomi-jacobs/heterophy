process IQTREE_PHYLOGENY {
    tag "${strategy}"
    label 'process_high'
    label 'process_long'

    conda "envs/phylo.yml"
    container "quay.io/biocontainers/iqtree:2.2.6--h21ec9f0_0"

    publishDir "${params.outdir}/trees/${strategy}", mode: 'copy'

    input:
    tuple val(strategy), path(alignment)

    output:
    tuple val(strategy), path("${strategy}.treefile"),          emit: tree
    tuple val(strategy), path("${strategy}.iqtree"),            emit: log
    tuple val(strategy), path("${strategy}.contree"),           emit: contree
    tuple val(strategy), path("${strategy}.ufboot"),            emit: bootstrap

    script:
    """
    iqtree2 \\
        -s ${alignment} \\
        -m ${params.iqtree_model} \\
        -B ${params.iqtree_bootstrap} \\
        -T ${params.iqtree_threads} \\
        --prefix ${strategy} \\
        --redo \\
        -nt AUTO
    """
}

/*
========================================================================================
    MODULE: tree_comparison — Compare topologies & bootstrap support between strategies
========================================================================================
*/

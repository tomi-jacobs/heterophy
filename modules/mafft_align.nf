process MAFFT_ALIGN {
    tag "${strategy}"
    label 'process_high'
    label 'process_long'

    conda "envs/phylo.yml"
    container "quay.io/biocontainers/mafft:7.520--hec16e2b_1"

    publishDir "${params.outdir}/alignments/mafft/${strategy}", mode: 'copy'

    input:
    tuple val(strategy), path(fastas)

    output:
    tuple val(strategy), path("${strategy}.aligned.fa"), emit: alignment

    script:
    """
    # Concatenate all species FASTAs
    cat ${fastas} > ${strategy}.combined.fa

    # Align
    mafft \\
        ${params.mafft_args} \\
        --thread ${task.cpus} \\
        ${strategy}.combined.fa > ${strategy}.aligned.fa
    """
}

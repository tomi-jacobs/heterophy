process RETAIN_HETEROZYGOSITY {
    tag "${meta.id}"
    label 'process_low'

    conda "envs/python.yml"
    container "heterophy/python:1.0.0"

    publishDir "${params.outdir}/strategies/retained", mode: 'copy'

    input:
    tuple val(meta), path(masked_fasta), path(filtered_vcf)

    output:
    tuple val(meta), path("${meta.id}.retained.fa"), emit: fasta

    script:
    """
    apply_heterozygosity_strategy.py \\
        --fasta      ${masked_fasta} \\
        --vcf        ${filtered_vcf} \\
        --strategy   retain \\
        --sample     ${meta.id} \\
        --output     ${meta.id}.retained.fa
    """
}

/*
========================================================================================
    MODULE: mafft_align — Multiple sequence alignment across all species
========================================================================================
*/

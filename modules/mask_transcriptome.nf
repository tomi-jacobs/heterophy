process MASK_TRANSCRIPTOME {
    tag "${meta.id}"
    label 'process_low'

    conda "envs/python.yml"
    container "heterophy/python:1.0.0"

    publishDir "${params.outdir}/masked_transcriptomes", mode: 'copy'

    input:
    tuple val(meta), path(transcriptome), path(coverage_table)

    output:
    tuple val(meta), path("${meta.id}.masked.fa"), emit: masked_fasta
    tuple val(meta), path("${meta.id}.masking_stats.json"), emit: masking_stats

    script:
    """
    mask_transcriptome.py \\
        --fasta       ${transcriptome} \\
        --coverage    ${coverage_table} \\
        --min_depth   ${params.min_coverage} \\
        --min_length  ${params.min_length} \\
        --sample      ${meta.id} \\
        --output      ${meta.id}.masked.fa \\
        --stats       ${meta.id}.masking_stats.json
    """
}

/*
========================================================================================
    MODULE: snp_calling — Call SNPs using BCFtools mpileup + call
========================================================================================
*/

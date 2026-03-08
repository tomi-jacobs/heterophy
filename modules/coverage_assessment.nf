process COVERAGE_ASSESSMENT {
    tag "${meta.id}"
    label 'process_medium'

    conda "envs/alignment.yml"
    container "quay.io/biocontainers/samtools:1.18--h50ea8bc_1"

    publishDir "${params.outdir}/coverage/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai), path(transcriptome)

    output:
    tuple val(meta), path("${meta.id}.coverage.tsv"),   emit: coverage_table
    tuple val(meta), path("${meta.id}.low_cov.bed"),    emit: low_cov_bed
    tuple val(meta), path("${meta.id}.coverage_stats.json"), emit: stats

    script:
    """
    # Per-nucleotide depth
    samtools depth \\
        -a \\
        -Q ${params.min_map_quality} \\
        ${bam} > ${meta.id}.coverage.tsv

    # Identify low-coverage regions
    coverage_assessment.py \\
        --coverage    ${meta.id}.coverage.tsv \\
        --min_depth   ${params.min_coverage} \\
        --sample      ${meta.id} \\
        --bed_out      ${meta.id}.low_cov.bed \\
        --stats_out    ${meta.id}.coverage_stats.json
    """
}

/*
========================================================================================
    MODULE: mask_transcriptome — Replace low-coverage positions with N
    RATIONALE: Unknown (N) states are treated as all-states in phylogenetic
               likelihood calculations — safer than retaining incorrect bases.
========================================================================================
*/

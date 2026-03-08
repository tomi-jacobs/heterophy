process COLLAPSE_HETEROZYGOSITY {
    tag "${meta.id}"
    label 'process_low'

    conda "envs/python.yml"
    container "heterophy/python:1.0.0"

    publishDir "${params.outdir}/strategies/collapsed", mode: 'copy'

    input:
    tuple val(meta), path(masked_fasta), path(filtered_vcf)

    output:
    tuple val(meta), path("${meta.id}.collapsed.fa"), emit: fasta

    script:
    """
    apply_heterozygosity_strategy.py \\
        --fasta      ${masked_fasta} \\
        --vcf        ${filtered_vcf} \\
        --strategy   collapse \\
        --sample     ${meta.id} \\
        --output     ${meta.id}.collapsed.fa
    """
}

/*
========================================================================================
    MODULE: retain_heterozygosity — Strategy 2: IUPAC ambiguity codes at het sites
    RATIONALE: Preserves allelic variation signaling recent divergence.
               Risk: may introduce non-existent genotypes into downstream inference.
========================================================================================
*/

process SNP_FILTERING {
    tag "${meta.id}"
    label 'process_low'

    conda "envs/variants.yml"
    container "quay.io/biocontainers/bcftools:1.18--h8b25389_0"

    publishDir "${params.outdir}/vcf/filtered/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(raw_vcf)

    output:
    tuple val(meta), path("${meta.id}.filtered.vcf.gz"), emit: filtered_vcf
    tuple val(meta), path("${meta.id}.snp_stats.json"),  emit: stats

    script:
    """
    # Quality filtering
    bcftools filter \\
        --include 'QUAL>=30 && FORMAT/DP>=${params.min_snp_depth} && MQ>=${params.min_map_quality}' \\
        --output-type z \\
        --output ${meta.id}.filtered.vcf.gz \\
        ${raw_vcf}

    bcftools index ${meta.id}.filtered.vcf.gz

    # Generate SNP stats for report
    snp_stats.py \\
        --vcf    ${meta.id}.filtered.vcf.gz \\
        --sample ${meta.id} \\
        --output ${meta.id}.snp_stats.json
    """
}

/*
========================================================================================
    MODULE: collapse_heterozygosity — Strategy 1: consensus base at het sites
    RATIONALE: Simplifies allelic variation for traditional phylogenetic models.
               Information is lost but avoids introduction of non-existent genotypes.
========================================================================================
*/

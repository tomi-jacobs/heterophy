process SNP_CALLING {
    tag "${meta.id}"
    label 'process_medium'

    conda "envs/variants.yml"
    container "quay.io/biocontainers/bcftools:1.18--h8b25389_0"

    publishDir "${params.outdir}/vcf/raw/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(reads), path(masked_fasta)

    output:
    tuple val(meta), path("${meta.id}.raw.vcf.gz"), emit: raw_vcf

    script:
    def r1 = reads[0]
    def r2 = reads[1]
    """
    # Index masked FASTA for samtools/bcftools
    samtools faidx ${masked_fasta}

    # Build HISAT2 index from masked FASTA
    hisat2-build -p ${task.cpus} ${masked_fasta} masked_index

    # Align reads to masked transcriptome
    hisat2 \\
        -p ${task.cpus} \\
        -x masked_index \\
        -1 ${r1} \\
        -2 ${r2} \\
        --no-spliced-alignment \\
    | samtools sort -@ ${task.cpus} -o ${meta.id}.masked.sorted.bam

    samtools index ${meta.id}.masked.sorted.bam

    # SNP calling: mpileup → call
    bcftools mpileup \\
        --threads ${task.cpus} \\
        --fasta-ref ${masked_fasta} \\
        --min-BQ ${params.min_base_quality} \\
        --min-MQ ${params.min_map_quality} \\
        --annotate FORMAT/AD,FORMAT/DP \\
        ${meta.id}.masked.sorted.bam \\
    | bcftools call \\
        --multiallelic-caller \\
        --variants-only \\
        --output-type z \\
        --output ${meta.id}.raw.vcf.gz

    bcftools index ${meta.id}.raw.vcf.gz
    """
}

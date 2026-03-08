process SAMTOOLS_SORT_INDEX {
    tag "${meta.id}"
    label 'process_medium'

    conda "envs/alignment.yml"
    container "quay.io/biocontainers/samtools:1.18--h50ea8bc_1"

    publishDir "${params.outdir}/alignments/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.id}.sorted.bam"), path("${meta.id}.sorted.bam.bai"), emit: bam_bai

    script:
    """
    samtools sort \\
        -@ ${task.cpus} \\
        -o ${meta.id}.sorted.bam \\
        ${bam}

    samtools index ${meta.id}.sorted.bam
    """
}

/*
========================================================================================
    MODULE: coverage_assessment — Per-nucleotide coverage depth via Samtools
    OUTPUT: BED file of low-coverage regions + full coverage table
========================================================================================
*/

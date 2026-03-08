process HISAT2_ALIGN {
    tag "${meta.id}"
    label 'process_high'

    conda "envs/alignment.yml"
    container "quay.io/biocontainers/hisat2:2.2.1--h1b792b2_3"

    publishDir "${params.outdir}/alignments/${meta.id}", mode: 'copy', pattern: "*.log"

    input:
    tuple val(meta), path(reads), path(index_dir), path(transcriptome)

    output:
    tuple val(meta), path("${meta.id}.bam"),     emit: bam
    tuple val(meta), path("${meta.id}.hisat2.log"), emit: log

    script:
    def r1 = reads[0]
    def r2 = reads[1]
    """
    hisat2 \\
        -p ${task.cpus} \\
        -x ${index_dir}/${meta.id} \\
        -1 ${r1} \\
        -2 ${r2} \\
        --dta \\
        --no-spliced-alignment \\
        --summary-file ${meta.id}.hisat2.log \\
        2>> ${meta.id}.hisat2.log \\
    | samtools view -bS - > ${meta.id}.bam
    """
}

/*
========================================================================================
    MODULE: samtools_sort_index — Sort and index BAM file
========================================================================================
*/

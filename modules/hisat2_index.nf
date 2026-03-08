/*
========================================================================================
    MODULE: hisat2_index — Build splice-aware HISAT2 index from transcriptome FASTA
========================================================================================
*/
process HISAT2_INDEX {
    tag "${meta.id}"
    label 'process_medium'

    conda "envs/alignment.yml"
    container "quay.io/biocontainers/hisat2:2.2.1--h1b792b2_3"

    input:
    tuple val(meta), path(transcriptome)

    output:
    tuple val(meta), path("${meta.id}_index"), emit: index

    script:
    """
    mkdir -p ${meta.id}_index
    hisat2-build \\
        -p ${task.cpus} \\
        ${transcriptome} \\
        ${meta.id}_index/${meta.id}
    """
}

/*
========================================================================================
    MODULE: hisat2_align — Align raw reads back to transcriptome (splice-aware)
    RATIONALE: Splice-aware alignment is critical for eukaryotic transcriptomes.
               HISAT2 handles exon-exon junctions that Bowtie2 would mishandle.
========================================================================================
*/

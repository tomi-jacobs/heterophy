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
    # Index masked FASTA for BCFtools
    samtools faidx ${masked_fasta}

    # Re-align reads to masked transcriptome
    hisat2 \\
        -p ${task.cpus} \\
        -x ${masked_fasta} \\
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

/*
========================================================================================
    MODULE: snp_filtering — Filter SNPs for high-confidence calls
========================================================================================
*/
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
process MAFFT_ALIGN {
    tag "${strategy}"
    label 'process_high'

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

/*
========================================================================================
    MODULE: iqtree_phylogeny — Maximum likelihood phylogenetic inference
========================================================================================
*/
process IQTREE_PHYLOGENY {
    tag "${strategy}"
    label 'process_high'
    label 'process_long'

    conda "envs/phylo.yml"
    container "quay.io/biocontainers/iqtree:2.2.6--h21ec9f0_0"

    publishDir "${params.outdir}/trees/${strategy}", mode: 'copy'

    input:
    tuple val(strategy), path(alignment)

    output:
    tuple val(strategy), path("${strategy}.treefile"),          emit: tree
    tuple val(strategy), path("${strategy}.iqtree"),            emit: log
    tuple val(strategy), path("${strategy}.contree"),           emit: contree
    tuple val(strategy), path("${strategy}.ufboot"),            emit: bootstrap

    script:
    """
    iqtree2 \\
        -s ${alignment} \\
        -m ${params.iqtree_model} \\
        -B ${params.iqtree_bootstrap} \\
        -T ${params.iqtree_threads} \\
        --prefix ${strategy} \\
        --redo \\
        -nt AUTO
    """
}

/*
========================================================================================
    MODULE: tree_comparison — Compare topologies & bootstrap support between strategies
========================================================================================
*/
process TREE_COMPARISON {
    tag "Comparing trees"
    label 'process_low'

    conda "envs/python.yml"
    container "heterophy/python:1.0.0"

    publishDir "${params.outdir}/comparison", mode: 'copy'

    input:
    path trees   // All .treefile outputs collected

    output:
    path "comparison_results.json",           emit: comparison_json
    path "topology_comparison.pdf",           emit: topology_plot
    path "bootstrap_comparison.pdf",          emit: bootstrap_plot
    path "robinson_foulds_distances.tsv",     emit: rf_distances
    path "bootstrap_support_summary.tsv",     emit: support_summary

    script:
    """
    compare_trees.py \\
        --trees      ${trees} \\
        --outdir     . \\
        --output     comparison_results.json
    """
}

/*
========================================================================================
    MODULE: generate_report — HTML report summarising the full analysis
========================================================================================
*/
process GENERATE_REPORT {
    tag "Generating report"
    label 'process_low'

    conda "envs/python.yml"
    container "heterophy/python:1.0.0"

    publishDir "${params.outdir}/report", mode: 'copy'

    input:
    path inputs   // coverage tables, SNP stats, comparison JSON

    output:
    path "heterophy_report.html", emit: report

    script:
    """
    generate_report.py \\
        --outdir       . \\
        --title        "${params.report_title}" \\
        --pipeline_ver "1.0.0"
    """
}

#!/usr/bin/env nextflow

/*
========================================================================================
    HeteroPhy: SNP & Heterozygosity-Aware Phylogenetic Inference Pipeline
========================================================================================
    GitHub    : https://github.com/[your-lab]/heterophy
    Author    : [Your Name]
    Lab       : [Your Lab / University]
    License   : MIT
    Citation  : [Your paper citation once published]
----------------------------------------------------------------------------------------
    A Nextflow DSL2 pipeline for investigating the role of SNPs and heterozygosity
    in resolving evolutionary relationships from de novo transcriptome assemblies.
    Designed for non-model organisms with rapid radiation dynamics.
========================================================================================
*/

nextflow.enable.dsl = 2

// ─────────────────────────────────────────────────────────────────────────────
// PARAMETER DEFAULTS
// ─────────────────────────────────────────────────────────────────────────────



// ─────────────────────────────────────────────────────────────────────────────
// HELP MESSAGE
// ─────────────────────────────────────────────────────────────────────────────

def helpMessage() {
    log.info """
    ╔══════════════════════════════════════════════════════════════════╗
    ║              HeteroPhy v1.0.0 — Usage Information               ║
    ╚══════════════════════════════════════════════════════════════════╝

    USAGE:
        nextflow run main.nf --input samplesheet.csv --outdir results

    REQUIRED ARGUMENTS:
        --input         Path to samplesheet CSV (sample,reads_r1,reads_r2,transcriptome)
        --outdir        Output directory [default: results]

    COVERAGE & MASKING:
        --min_coverage      Minimum depth to retain nucleotide [default: 5]
        --min_length        Minimum post-masking contig length [default: 200]

    SNP CALLING:
        --min_base_quality  Minimum base quality score [default: 20]
        --min_map_quality   Minimum mapping quality [default: 20]
        --min_snp_depth     Minimum read depth for SNP call [default: 10]
        --max_missing       Maximum missing data fraction per site [default: 0.5]

    PHYLOGENETICS:
        --mafft_args        MAFFT alignment arguments [default: --auto]
        --iqtree_model      IQ-TREE substitution model [default: TEST]
        --iqtree_bootstrap  Bootstrap replicates [default: 1000]

    PROFILES:
        -profile conda          Use Conda environments
        -profile docker         Use Docker containers
        -profile singularity    Use Singularity containers
        -profile slurm          Run on SLURM HPC cluster
        -profile cloud          Run on AWS/GCP

    EXAMPLES:
        # Local run with conda
        nextflow run main.nf --input samplesheet.csv -profile conda

        # HPC cluster with singularity
        nextflow run main.nf --input samplesheet.csv -profile slurm,singularity

        # Cloud run
        nextflow run main.nf --input samplesheet.csv -profile cloud,docker

    SAMPLESHEET FORMAT (CSV):
        sample,reads_r1,reads_r2,transcriptome
        species1,/path/to/sp1_R1.fastq.gz,/path/to/sp1_R2.fastq.gz,/path/to/sp1.fa
        species2,/path/to/sp2_R1.fastq.gz,/path/to/sp2_R2.fastq.gz,/path/to/sp2.fa
    """.stripIndent()
}

if (params.help) {
    helpMessage()
    exit 0
}

// ─────────────────────────────────────────────────────────────────────────────
// VALIDATE INPUTS
// ─────────────────────────────────────────────────────────────────────────────

if (!params.input) {
    error "ERROR: --input samplesheet is required. Run with --help for usage."
}

// ─────────────────────────────────────────────────────────────────────────────
// IMPORT MODULES
// ─────────────────────────────────────────────────────────────────────────────

include { VALIDATE_SAMPLESHEET     } from './modules/validate_samplesheet'
include { HISAT2_INDEX             } from './modules/hisat2_index'
include { HISAT2_ALIGN             } from './modules/hisat2_align'
include { SAMTOOLS_SORT_INDEX      } from './modules/samtools_sort_index'
include { COVERAGE_ASSESSMENT      } from './modules/coverage_assessment'
include { MASK_TRANSCRIPTOME       } from './modules/mask_transcriptome'
include { SNP_CALLING              } from './modules/snp_calling'
include { SNP_FILTERING            } from './modules/snp_filtering'
include { COLLAPSE_HETEROZYGOSITY  } from './modules/collapse_heterozygosity'
include { RETAIN_HETEROZYGOSITY    } from './modules/retain_heterozygosity'
include { MAFFT_ALIGN              } from './modules/mafft_align'
include { IQTREE_PHYLOGENY         } from './modules/iqtree_phylogeny'
include { TREE_COMPARISON          } from './modules/tree_comparison'
include { GENERATE_REPORT          } from './modules/generate_report'

// ─────────────────────────────────────────────────────────────────────────────
// MAIN WORKFLOW
// ─────────────────────────────────────────────────────────────────────────────

workflow {

    // ── Banner ────────────────────────────────────────────────────────────────
    log.info """
    ╔══════════════════════════════════════════════════════════════════╗
    ║    H E T E R O P H Y   v 1 . 0 . 0                              ║
    ║    SNP & Heterozygosity-Aware Phylogenetic Inference             ║
    ╚══════════════════════════════════════════════════════════════════╝
    input        : ${params.input}
    outdir       : ${params.outdir}
    min_coverage : ${params.min_coverage}
    iqtree_model : ${params.iqtree_model}
    bootstrap    : ${params.iqtree_bootstrap}
    ─────────────────────────────────────────────────────────────────
    """.stripIndent()

    // ── Step 0: Parse & validate samplesheet ─────────────────────────────────
    ch_samplesheet = Channel
        .fromPath(params.input, checkIfExists: true)

    VALIDATE_SAMPLESHEET(ch_samplesheet)

    ch_samples = VALIDATE_SAMPLESHEET.out.validated
        .splitCsv(header: true)
        .map { row ->
            def meta = [id: row.sample]
            def reads = [
                file(row.reads_r1, checkIfExists: true),
                file(row.reads_r2, checkIfExists: true)
            ]
            def transcriptome = file(row.transcriptome, checkIfExists: true)
            return [meta, reads, transcriptome]
        }

    ch_reads         = ch_samples.map { meta, reads, tx -> [meta, reads] }
    ch_transcriptomes = ch_samples.map { meta, reads, tx -> [meta, tx] }

    // ── Step 1: Build HISAT2 index for each transcriptome ────────────────────
    HISAT2_INDEX(ch_transcriptomes)

    // ── Step 2: Align reads back to transcriptome (splice-aware) ─────────────
    ch_align_input = ch_reads
        .join(HISAT2_INDEX.out.index, by: [0])
        .join(ch_transcriptomes, by: [0])

    HISAT2_ALIGN(ch_align_input)

    // ── Step 3: Sort & index BAM files ───────────────────────────────────────
    SAMTOOLS_SORT_INDEX(HISAT2_ALIGN.out.bam)

    // ── Step 4: Nucleotide-level coverage assessment ─────────────────────────
    ch_coverage_input = SAMTOOLS_SORT_INDEX.out.bam_bai
        .join(ch_transcriptomes, by: [0])

    COVERAGE_ASSESSMENT(ch_coverage_input)

    // ── Step 5: Mask low-coverage regions → masked FASTA ────────────────────
    ch_mask_input = ch_transcriptomes
        .join(COVERAGE_ASSESSMENT.out.coverage_table, by: [0])

    MASK_TRANSCRIPTOME(ch_mask_input)

    // ── Step 6: SNP calling on masked transcriptome ──────────────────────────
    ch_snp_input = ch_reads
        .join(MASK_TRANSCRIPTOME.out.masked_fasta, by: [0])

    SNP_CALLING(ch_snp_input)

    // ── Step 7: Filter SNPs for high confidence ──────────────────────────────
    SNP_FILTERING(SNP_CALLING.out.raw_vcf)

    // ── Step 8a: Collapse heterozygous sites ─────────────────────────────────
    ch_collapse_input = MASK_TRANSCRIPTOME.out.masked_fasta
        .join(SNP_FILTERING.out.filtered_vcf, by: [0])

    if (params.run_collapsed) {
        COLLAPSE_HETEROZYGOSITY(ch_collapse_input)
    }

    // ── Step 8b: Retain heterozygous sites (IUPAC ambiguity codes) ───────────
    if (params.run_retained) {
        RETAIN_HETEROZYGOSITY(ch_collapse_input)
    }

    // ── Step 9: Multiple sequence alignment with MAFFT ───────────────────────
    // Collect all species FASTAs per strategy into one channel for alignment
    if (params.run_collapsed && params.run_retained) {
        ch_collapsed_fasta = COLLAPSE_HETEROZYGOSITY.out.fasta
            .map { meta, fasta -> ["collapsed", meta.id, fasta] }
        ch_retained_fasta  = RETAIN_HETEROZYGOSITY.out.fasta
            .map { meta, fasta -> ["retained", meta.id, fasta] }

        ch_for_alignment = ch_collapsed_fasta.mix(ch_retained_fasta)
            .groupTuple(by: 0)
            .map { strategy, ids, fastas -> [strategy, fastas] }
    } else if (params.run_collapsed) {
        ch_for_alignment = COLLAPSE_HETEROZYGOSITY.out.fasta
            .collect()
            .map { fastas -> ["collapsed", fastas] }
    } else {
        ch_for_alignment = RETAIN_HETEROZYGOSITY.out.fasta
            .collect()
            .map { fastas -> ["retained", fastas] }
    }

    MAFFT_ALIGN(ch_for_alignment)

    // ── Step 10: Phylogenetic inference with IQ-TREE ─────────────────────────
    IQTREE_PHYLOGENY(MAFFT_ALIGN.out.alignment)

    // ── Step 11: Compare tree topologies & support values ────────────────────
    ch_trees_for_comparison = IQTREE_PHYLOGENY.out.tree
        .groupTuple(by: [])
        .collect()

    TREE_COMPARISON(IQTREE_PHYLOGENY.out.tree.collect())

    // ── Step 12: Generate HTML report ────────────────────────────────────────
    ch_report_inputs = COVERAGE_ASSESSMENT.out.coverage_table.collect()
        .combine(SNP_FILTERING.out.stats.collect())
        .combine(TREE_COMPARISON.out.comparison_json)

    GENERATE_REPORT(ch_report_inputs)

    // ── Completion message ────────────────────────────────────────────────────
    workflow.onComplete {
        log.info """
        ════════════════════════════════════════════════════════
        HeteroPhy pipeline complete!
        Status    : ${workflow.success ? 'SUCCESS ✓' : 'FAILED ✗'}
        Duration  : ${workflow.duration}
        Results   : ${params.outdir}
        Report    : ${params.outdir}/report/heterophy_report.html
        ════════════════════════════════════════════════════════
        """.stripIndent()
    }
}

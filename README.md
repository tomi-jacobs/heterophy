# HeteroPhy

**SNP & Heterozygosity-Aware Phylogenetic Inference Pipeline**

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A523.04.0-brightgreen.svg)](https://www.nextflow.io/)
[![run with conda](https://img.shields.io/badge/run%20with-conda-3EB049?logo=anaconda)](https://docs.conda.io/en/latest/miniconda.html)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg)](https://sylabs.io/docs/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

---

## Overview

**HeteroPhy** is a reproducible, scalable Nextflow pipeline that investigates the role of **Single Nucleotide Polymorphisms (SNPs)** and **heterozygosity** in resolving evolutionary relationships from de novo assembled transcriptomes.

The pipeline is designed for non-model organisms, particularly those undergoing rapid evolutionary radiation, where the phylogenetic signal resides in a small number of genetic differences and the influence of heterozygosity on inference has rarely been empirically tested.

### The core question HeteroPhy answers:

> **Does how we represent heterozygous sites in transcriptomic data change the evolutionary relationships we infer — and does it affect our confidence in those relationships?**

---

## Background

De novo transcriptome assembly of non-model organisms is rarely preceded by inbreeding, meaning assembled ORFs frequently contain a mixture of alleles. This raises a fundamental but largely unexplored question: does the choice of how to handle heterozygous sites — collapsing them to one allele versus retaining ambiguity — influence phylogenetic topology or statistical support?

HeteroPhy empirically tests this by:
1. Mapping reads back to assembled transcriptomes to calculate per-nucleotide coverage
2. Masking low-confidence positions to avoid biased inference
3. Calling SNPs on masked transcriptomes
4. Applying two strategies: **Collapse** (reference allele) and **Retain** (IUPAC codes)
5. Inferring independent phylogenies under each strategy
6. Quantifying the topological and support differences between them

---

## Pipeline Overview

```
Raw Reads + Assembled Transcriptomes
           │
           ▼
    ┌─────────────────────────────────┐
    │  Step 1: HISAT2 Index           │  Splice-aware index (eukaryotic)
    │  Step 2: HISAT2 Align           │  Read → transcriptome alignment
    │  Step 3: Samtools Sort/Index    │  BAM processing
    └────────────┬────────────────────┘
                 │
                 ▼
    ┌─────────────────────────────────┐
    │  Step 4: Coverage Assessment    │  Per-nucleotide depth (samtools depth)
    │  Step 5: Mask Transcriptome     │  Low-depth positions → N
    └────────────┬────────────────────┘
                 │
                 ▼
    ┌─────────────────────────────────┐
    │  Step 6: SNP Calling            │  BCFtools mpileup + call
    │  Step 7: SNP Filtering          │  QUAL≥30, DP≥10, MQ≥20
    └────────────┬────────────────────┘
                 │
           ┌─────┴──────┐
           ▼             ▼
    ┌──────────┐   ┌───────────┐
    │ COLLAPSE │   │  RETAIN   │
    │ Strategy │   │ Strategy  │
    │ Ref allele│  │IUPAC codes│
    └─────┬────┘   └─────┬─────┘
          │               │
          ▼               ▼
    ┌─────────────────────────────────┐
    │  Step 9:  MAFFT Alignment       │  Per strategy
    │  Step 10: IQ-TREE2 Phylogeny    │  ML + UFBoot
    └────────────┬────────────────────┘
                 │
                 ▼
    ┌─────────────────────────────────┐
    │  Step 11: Tree Comparison       │  RF distance, KS test
    │  Step 12: HTML Report           │  Publication-ready
    └─────────────────────────────────┘
```

---

## Key Design Decisions

### Why HISAT2 (not Bowtie2)?

HISAT2 is splice-aware, essential for eukaryotic RNA-seq data. Eukaryotic transcriptomes contain exon-exon junctions; Bowtie2 would misalign reads spanning these junctions, artificially reducing coverage at junction sites and misclassifying them as low-coverage.

### Why mask to N (not remove)?

In phylogenetic likelihood models, an `N` (unknown state) is treated as all states simultaneously — it contributes no directional information and thus introduces no bias. A nucleotide with incorrect identity, by contrast, actively misleads topology inference. This is especially critical for rapidly radiating species where the phylogenetic signal is a small number of genuine differences between closely related species.

### Collapse vs Retain — what do each mean?

| Strategy | Heterozygous site treatment | Risk | Benefit |
|---|---|---|---|
| **Collapse** | Use reference allele | Information loss | Compatible with standard substitution models |
| **Retain** | IUPAC ambiguity code | Phantom genotypes | Preserves recent divergence signal |

---

## Installation

### Prerequisites

- [Nextflow](https://www.nextflow.io/) ≥ 23.04.0
- One of: [Conda/Mamba](https://docs.conda.io), [Docker](https://docker.com), or [Singularity](https://sylabs.io)

### Install Nextflow

```bash
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
```

### Clone the repository

```bash
git clone https://github.com/tomi-jacobs/heterophy.git
cd heterophy
```

---

## Quick Start

### 1. Prepare your samplesheet

Create a CSV file with the following columns:

```csv
sample,reads_r1,reads_r2,transcriptome
species1,/path/to/sp1_R1.fastq.gz,/path/to/sp1_R2.fastq.gz,/path/to/sp1_cds.fa
species2,/path/to/sp2_R1.fastq.gz,/path/to/sp2_R2.fastq.gz,/path/to/sp2_cds.fa
```

> **Note:** Transcriptome FASTAs should contain CDS sequences (coding sequences), ideally from a Semblans or Trinity assembly. Raw transcripts are also supported.

### 2. Run the pipeline

```bash
# Local run with conda
nextflow run main.nf \
    --input samplesheet.csv \
    --outdir results \
    -profile conda

# HPC cluster with Singularity
nextflow run main.nf \
    --input samplesheet.csv \
    --outdir results \
    -profile slurm,singularity

# Test with provided test data
nextflow run main.nf \
    -profile test,conda
```

---

## Parameters

### Required

| Parameter | Description |
|---|---|
| `--input` | Path to samplesheet CSV |
| `--outdir` | Output directory |

### Coverage & Masking

| Parameter | Default | Description |
|---|---|---|
| `--min_coverage` | `5` | Minimum read depth to retain a nucleotide position |
| `--min_length` | `200` | Minimum retained contig length post-masking |

### SNP Calling

| Parameter | Default | Description |
|---|---|---|
| `--min_base_quality` | `20` | Minimum base quality for SNP calling |
| `--min_map_quality` | `20` | Minimum mapping quality |
| `--min_snp_depth` | `10` | Minimum depth for a SNP call |

### Phylogenetics

| Parameter | Default | Description |
|---|---|---|
| `--mafft_args` | `--auto` | MAFFT strategy |
| `--iqtree_model` | `TEST` | ModelFinder model selection |
| `--iqtree_bootstrap` | `1000` | Ultrafast bootstrap replicates |

---

## Output Structure

```
results/
├── alignments/
│   ├── {sample}/           # HISAT2 BAM files and logs
│   └── mafft/              # Multiple sequence alignments
│       ├── collapsed/
│       └── retained/
├── coverage/
│   └── {sample}/           # Coverage TSVs, BED files, stats JSON
├── masked_transcriptomes/  # Per-sample masked FASTAs
├── vcf/
│   ├── raw/{sample}/       # Pre-filter VCF
│   └── filtered/{sample}/  # High-confidence SNP VCF
├── strategies/
│   ├── collapsed/          # Collapsed heterozygosity FASTAs
│   └── retained/           # IUPAC-encoded FASTAs
├── trees/
│   ├── collapsed/          # IQ-TREE output (collapsed strategy)
│   └── retained/           # IQ-TREE output (retained strategy)
├── comparison/
│   ├── comparison_results.json
│   ├── robinson_foulds_distances.tsv
│   ├── bootstrap_support_summary.tsv
│   ├── topology_comparison.pdf
│   └── bootstrap_comparison.pdf
├── report/
│   └── heterophy_report.html   # ← Main output — open this
└── pipeline_info/
    ├── execution_report.html
    ├── execution_timeline.html
    └── execution_trace.txt
```

---

## Interpreting Results

### Robinson-Foulds Distance

| Normalized RF | Interpretation |
|---|---|
| 0.00 | Identical topologies — strategy has no effect |
| < 0.05 | Near-identical — strategy has negligible effect |
| 0.05 – 0.20 | Minor differences — strategy has modest influence |
| 0.20 – 0.50 | Moderate differences — strategy influences inference |
| > 0.50 | Major differences — strategy substantially alters phylogenetic inference |

### Bootstrap KS Test

A significant KS test (p < 0.05) indicates that the two strategies produce meaningfully different confidence distributions across nodes — one strategy may produce systematically better- or worse-supported trees.

---

## Citation

If you use HeteroPhy in your research, please cite:

```
PENDING!
[Tomi Jacobs] et al. (2025). HeteroPhy: A bioinformatics pipeline for investigating
the role of SNPs and heterozygosity in resolving evolutionary relationships from
de novo transcriptomes. [Journal]. doi: [pending]
```

### Tool citations (please also cite)

- **HISAT2**: Kim D et al. (2019) Graph-based genome alignment and genotyping with HISAT2. *Nature Methods*
- **Samtools/BCFtools**: Danecek P et al. (2021) Twelve years of SAMtools and BCFtools. *GigaScience*
- **MAFFT**: Katoh K & Standley DM (2013) MAFFT Multiple Sequence Alignment Software. *Molecular Biology and Evolution*
- **IQ-TREE2**: Minh BQ et al. (2020) IQ-TREE 2: New Models and Methods for Phylogenetic Inference. *Molecular Biology and Evolution*
- **ETE3**: Huerta-Cepas J et al. (2016) ETE 3: Reconstruction, Analysis, and Visualization of Phylogenomic Data. *Molecular Biology and Evolution*

---

## Contributing

We welcome contributions, bug reports, and feature requests via [GitHub Issues](https://github.com/[your-lab]/heterophy/issues).

---

## License

MIT License — see [LICENSE](LICENSE) for details.

---

## Contact

[Tomi Jacobs] · [Walker Lab] · [University of Illinois, Chicago]  
[tomijacobs.e@gmail.edu]

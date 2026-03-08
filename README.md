# HeteroPhy

**SNP & Heterozygosity-Aware Phylogenetic Inference Pipeline**

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A523.04.0-brightgreen.svg)](https://www.nextflow.io/)
[![run with conda](https://img.shields.io/badge/run%20with-conda-3EB049?logo=anaconda)](https://docs.conda.io/en/latest/miniconda.html)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg)](https://sylabs.io/docs/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

**Author:** Tomi Jacobs · Walker Lab  
**Repository:** https://github.com/tomi-jacobs/heterophy  
**Citation:** Jacobs T. et al. (in prep). HeteroPhy: A pipeline for investigating the role of SNPs and heterozygosity in resolving evolutionary relationships from de novo transcriptomes.

---

## The Problem HeteroPhy Solves

When working with non-model organisms — plants, animals, fungi — de novo transcriptome assembly is rarely preceded by inbreeding. This means assembled transcripts frequently contain a mixture of alleles from both copies of each gene. Every time a researcher builds a phylogeny from such data, they face a silent, largely untested decision: **what do you do with heterozygous sites?**

The two dominant approaches in practice are:

- **Collapse** — resolve heterozygous positions to a single allele (typically the reference), discarding the allelic variation entirely
- **Retain** — encode heterozygous positions using IUPAC ambiguity codes (e.g. R = A/G, Y = C/T), preserving both alleles in the alignment

Both approaches are used in the literature. Neither has been systematically tested to determine whether the choice changes the topology of the inferred phylogeny or the statistical confidence placed in its nodes. **For organisms undergoing rapid evolutionary radiation — where the phylogenetic signal consists of a small number of genuine differences between closely related species — this question is not academic. It is central to the reliability of every downstream conclusion.**

HeteroPhy is the first pipeline built specifically to answer it.

---

## Pipeline Overview

```
Raw Reads + De Novo Assembled Transcriptomes (any eukaryote)
                        │
                        ▼
    ┌───────────────────────────────────────┐
    │  Step 1: HISAT2 Index                 │  Splice-aware index
    │  Step 2: HISAT2 Align                 │  Reads → transcriptome
    │  Step 3: Samtools Sort & Index        │  BAM processing
    └──────────────────┬────────────────────┘
                       │
                       ▼
    ┌───────────────────────────────────────┐
    │  Step 4: Coverage Assessment          │  Per-nucleotide depth
    │  Step 5: Mask Transcriptome           │  Low-depth positions → N
    └──────────────────┬────────────────────┘
                       │
                       ▼
    ┌───────────────────────────────────────┐
    │  Step 6: SNP Calling                  │  BCFtools mpileup + call
    │  Step 7: SNP Filtering                │  QUAL≥30, DP≥10, MQ≥20
    └──────────────────┬────────────────────┘
                       │
              ┌────────┴────────┐
              ▼                 ▼
       ┌────────────┐    ┌────────────┐
       │  COLLAPSE  │    │   RETAIN   │
       │ Ref allele │    │IUPAC codes │
       └─────┬──────┘    └─────┬──────┘
             │                 │
             ▼                 ▼
    ┌───────────────────────────────────────┐
    │  Step 9:  MAFFT Alignment             │  Per strategy
    │  Step 10: IQ-TREE2 Phylogeny          │  ML + UFBoot
    └──────────────────┬────────────────────┘
                       │
                       ▼
    ┌───────────────────────────────────────┐
    │  Step 11: Tree Comparison             │  RF distance, KS test
    │  Step 12: HTML Report                 │  Interactive summary
    └───────────────────────────────────────┘
```

### Why HISAT2 (not Bowtie2)?

HISAT2 is splice-aware. Eukaryotic transcriptomes contain exon-exon junctions; reads spanning these junctions would be misaligned by Bowtie2, artificially reducing coverage at junction sites and causing them to be incorrectly masked. This distinction matters especially for organisms with complex gene structures.

### Why mask to N (not remove)?

In phylogenetic likelihood models, `N` (unknown state) is treated as all states simultaneously — it contributes no directional information and introduces no bias. A nucleotide with an incorrect identity, by contrast, actively misleads topology inference. Masking rather than removing preserves alignment length while eliminating the risk of false signal. This is especially critical for rapidly radiating species where the phylogenetic signal is a small number of genuine differences.

### Collapse vs. Retain — what each strategy means

| Strategy | Het site treatment | Risk | Benefit |
|---|---|---|---|
| **Collapse** | Reference allele used | Information loss | Compatible with standard substitution models |
| **Retain** | IUPAC ambiguity code | Phantom genotypes possible | Preserves recent divergence signal |

---

## Installation

### Prerequisites

- [Nextflow](https://www.nextflow.io/) ≥ 23.04.0
- One of: [Conda/Mamba](https://docs.conda.io), [Docker](https://docker.com), or [Singularity](https://sylabs.io)

### Install Nextflow

```bash
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
nextflow -version
```

### Clone the repository

```bash
git clone https://github.com/tomi-jacobs/heterophy.git
cd heterophy
```

---

## Quick Start

### 1. Prepare your samplesheet

```csv
sample,reads_r1,reads_r2,transcriptome
species1,/path/to/sp1_R1.fastq.gz,/path/to/sp1_R2.fastq.gz,/path/to/sp1.cds.fa
species2,/path/to/sp2_R1.fastq.gz,/path/to/sp2_R2.fastq.gz,/path/to/sp2.cds.fa
```

Transcriptome FASTAs should be CDS sequences from a de novo assembler (e.g. Trinity + TransDecoder, or Semblans). Raw transcripts are also supported. Works with any eukaryote — plants, animals, fungi.

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

# Test run on 3 species with reduced bootstraps
nextflow run main.nf \
    --input samplesheet_small.csv \
    --outdir results_test \
    --iqtree_bootstrap 100 \
    -profile conda
```

Always use `-resume` to restart from the last completed step if anything fails:

```bash
nextflow run main.nf --input samplesheet.csv --outdir results -profile conda -resume
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
| `--min_length` | `200` | Minimum contig length post-masking |

### SNP Calling

| Parameter | Default | Description |
|---|---|---|
| `--min_base_quality` | `20` | Minimum base quality for SNP calling |
| `--min_map_quality` | `20` | Minimum mapping quality |
| `--min_snp_depth` | `10` | Minimum depth for a SNP call |

### Phylogenetics

| Parameter | Default | Description |
|---|---|---|
| `--mafft_args` | `--auto` | MAFFT alignment strategy |
| `--iqtree_model` | `TEST` | ModelFinder model selection |
| `--iqtree_bootstrap` | `1000` | Ultrafast bootstrap replicates |

### Execution

| Parameter | Default | Description |
|---|---|---|
| `--max_cpus` | `16` | Maximum CPUs per process |
| `--max_memory` | `128.GB` | Maximum memory per process |
| `--max_time` | `480.h` | Maximum walltime per process |

---

## Output Structure

```
results/
├── alignments/
│   ├── {sample}/           # HISAT2 BAM files and alignment logs
│   └── mafft/              # Multiple sequence alignments
│       ├── collapsed/
│       └── retained/
├── coverage/
│   └── {sample}/           # Per-nucleotide coverage TSV, BED, stats JSON
├── masked_transcriptomes/  # Per-sample masked FASTAs
├── vcf/
│   ├── raw/{sample}/       # Pre-filter VCF
│   └── filtered/{sample}/  # High-confidence SNP VCF + stats
├── strategies/
│   ├── collapsed/          # Collapsed heterozygosity FASTAs
│   └── retained/           # IUPAC-encoded FASTAs
├── trees/
│   ├── collapsed/          # IQ-TREE2 output — collapsed strategy
│   └── retained/           # IQ-TREE2 output — retained strategy
├── comparison/
│   ├── comparison_results.json
│   ├── robinson_foulds_distances.tsv
│   ├── bootstrap_support_summary.tsv
│   ├── topology_comparison.pdf
│   └── bootstrap_comparison.pdf
├── report/
│   └── heterophy_report.html   # ← Start here
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
| 0.00 | Identical — heterozygosity treatment has no effect on topology |
| < 0.05 | Near-identical — negligible topological effect |
| 0.05–0.20 | Minor differences — modest influence |
| 0.20–0.50 | Moderate differences — heterozygosity influences relationship inference |
| > 0.50 | Major differences — heterozygosity substantially alters inferred phylogeny |

### Bootstrap KS Test

A significant KS test (p < 0.05) indicates that the two strategies produce meaningfully different confidence distributions across nodes — one strategy yields systematically better- or worse-supported trees for your dataset.

### Recommendation

If RF = 0 and KS is non-significant: both strategies yield equivalent results; report either.

If RF > 0 or KS is significant: report both trees, discuss discordant nodes explicitly, and treat the comparison itself as a finding — your dataset is sensitive to heterozygosity treatment.

---

## Execution Profiles

| Profile | Use case |
|---|---|
| `conda` | Local run using Conda/Mamba environments |
| `docker` | Containerised run with Docker |
| `singularity` | Containerised run with Singularity (HPC) |
| `slurm` | SLURM HPC cluster |
| `pbs` | PBS/Torque HPC cluster |
| `aws` | AWS Batch cloud execution |
| `gcp` | Google Cloud Life Sciences |

Profiles can be combined: `-profile slurm,singularity`

---

## Citation

If you use HeteroPhy in your research, please cite:

```
<<<<<<< HEAD
PENDING!
[Tomi Jacobs] et al. (2025). HeteroPhy: A bioinformatics pipeline for investigating
the role of SNPs and heterozygosity in resolving evolutionary relationships from
de novo transcriptomes. [Journal]. doi: [pending]
=======
Jacobs T. et al. (in prep). HeteroPhy: A pipeline for investigating the role of SNPs
and heterozygosity in resolving evolutionary relationships from de novo transcriptomes.
Walker Lab.
>>>>>>> 1736a16 (Rewrite README: clean problem framing, remove competitor comparisons)
```

### Tool citations

Please also cite the tools HeteroPhy depends on:

- **HISAT2**: Kim D et al. (2019) Graph-based genome alignment and genotyping with HISAT2. *Nature Methods* 16:189–191
- **Samtools/BCFtools**: Danecek P et al. (2021) Twelve years of SAMtools and BCFtools. *GigaScience* 10:giab008
- **MAFFT**: Katoh K & Standley DM (2013) MAFFT Multiple Sequence Alignment Software. *Molecular Biology and Evolution* 30:772–780
- **IQ-TREE2**: Minh BQ et al. (2020) IQ-TREE 2: New Models and Methods for Phylogenetic Inference. *Molecular Biology and Evolution* 37:1530–1534
- **ETE3**: Huerta-Cepas J et al. (2016) ETE 3: Reconstruction, Analysis, and Visualization of Phylogenomic Data. *Molecular Biology and Evolution* 33:1635–1638
- **Nextflow**: Di Tommaso P et al. (2017) Nextflow enables reproducible computational workflows. *Nature Biotechnology* 35:316–319

---

## Contributing

Bug reports, feature requests, and pull requests are welcome via [GitHub Issues](https://github.com/tomi-jacobs/heterophy/issues).

---

## License

MIT License — see [LICENSE](LICENSE) for details.

---

## Contact

<<<<<<< HEAD
[Tomi Jacobs] · [Walker Lab] · [University of Illinois, Chicago]  
[tomijacobs.e@gmail.edu]
=======
Tomi Jacobs · Walker Lab  
GitHub: [@tomi-jacobs](https://github.com/tomi-jacobs)
>>>>>>> 1736a16 (Rewrite README: clean problem framing, remove competitor comparisons)

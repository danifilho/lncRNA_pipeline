# IWGC lncRNA Prediction Pipeline

A reproducible **Snakemake + Singularity/Apptainer** workflow for *de novo* discovery of
long non-coding RNAs (lncRNAs) in non-model plants, expanding existing genome annotations.

The pipeline assembles a reference-guided transcriptome from public RNA-seq data, extracts
novel transcript candidates, and keeps only those that **four independent coding-potential
tools agree are non-coding** and that have **no protein hit** — a high-confidence consensus
set of lncRNAs.

> Adapted from the Plant-LncRNA-pipeline-v2 by Xue-Chan Tian *et al.*
> ([DOI](https://doi.org/10.1093/hr/uhae041) ·
> [GitHub](https://github.com/xuechantian/Plant-LncRNA-pipeline-v2)) and
> [FEELnc](https://github.com/tderrien/FEELnc), re-implemented as a containerized Snakemake
> workflow for HPC (tested on the MSU HPCC / ICER).

---

## Table of contents

- [Overview](#overview)
- [Workflow diagram](#workflow-diagram)
- [Inputs](#inputs)
- [Pipeline steps](#pipeline-steps)
- [Containers](#containers)
- [Outputs](#outputs)
- [How to run](#how-to-run)
- [Test results (Arabidopsis TAIR10)](#test-results-arabidopsis-tair10)
- [Repository layout](#repository-layout)
- [Further documentation](#further-documentation)

---

## Overview

| | |
|---|---|
| **Goal** | Discover high-confidence lncRNAs not present in a species' current annotation |
| **Engine** | [Snakemake](https://snakemake.github.io/) 9.x |
| **Reproducibility** | One [Singularity/Apptainer](https://apptainer.org/) image per tool (15 images) |
| **Inputs** | Reference genome + annotation, RNA-seq reads, a protein database |
| **Output** | A consensus set of lncRNA IDs, GTF, FASTA, and FEELnc classifications |
| **Consensus rule** | `FEELnc ∩ CPAT ∩ LncFinder ∩ PlantLncBoost ∩ (no DIAMOND protein hit)` |

The workflow is fully described in [`workflow/Snakefile`](workflow/Snakefile); all paths and
parameters are configured in [`config/config.yaml`](config/config.yaml).

---

## Workflow diagram

```mermaid
flowchart TD
    subgraph IN[1 . Inputs]
        G[Reference genome FASTA]
        A[Annotation GTF]
        R[RNA-seq runs / FASTQ]
        U[UniProt protein DB]
    end

    R -->|prefetch / fasterq-dump| RAW[raw FASTQ]
    RAW -->|fastp| TRIM[trimmed reads]
    G --> IDX
    A --> IDX[STAR genome index]
    TRIM -->|STAR| BAM[sorted BAM]
    IDX --> BAM
    BAM -->|samtools index| BAI[indexed BAM]
    BAI -->|StringTie assemble| ASM[per-sample GTF]
    A --> MRG
    ASM -->|StringTie merge| MRG[merged transcriptome GTF]
    MRG -->|grep MSTRG + gffread| CAND[candidate GTF / FASTA]

    CAND --> FE[FEELnc_filter]
    CAND --> CP[CPAT]
    CAND --> LF[LncFinder]
    CAND --> LB[PlantLncBoost]
    CAND --> DI[DIAMOND blastx]
    U --> DI

    FE --> INT
    CP --> INT
    LF --> INT
    LB --> INT
    DI --> INT[Intersect IDs]
    INT --> FINAL[Consensus lncRNAs<br/>IDs · GTF · FASTA]
    FINAL -->|FEELnc_classifier| CLS[lncRNA classes]
    FINAL -.optional.->|minimap2 vs truth set| VAL[validation metrics]
```

The Snakefile expands to **22 jobs** for a 2-sample run (one `fastp`/`star_align`/
`samtools_index`/`stringtie_assemble` per sample, plus the shared downstream steps).

---

## Inputs

All inputs are declared in `config/config.yaml` (paths support `{species}`/`{prefix}`
templating).

| Input | Config key | Description | How to obtain |
|-------|-----------|-------------|---------------|
| **Reference genome** | `reference_fasta` | Cleaned chromosome FASTA (`{species}.chromosomes.fa`) | NCBI `datasets` → `scripts/ncbi_datasets_cleanup.py` |
| **Annotation** | `annotation_gtf` | Protein-coding annotation GTF (`{species}_mRNA.gtf`) | same cleanup script |
| **RNA-seq** | `sra_accessions` or `raw_fastq_dir` | SRA accessions to download, **or** pre-downloaded FASTQ | PLncDB / SRA / ENA |
| **Protein DB** | `diamond_db` / `uniprot_fasta` | DIAMOND database, or a SwissProt FASTA to build one | [UniProt](https://www.uniprot.org/) |
| **Truth set** *(optional)* | `truth_input_fasta` | Known lncRNAs for benchmarking (minimap2 step) | [PLncDB](https://www.tobaccodb.org/plncdb/) / EVlncRNA |

**Preparing the reference** (NCBI `datasets`):

```bash
datasets download genome taxon "Arabidopsis thaliana" --reference --include genome,gff3
unzip ncbi_dataset.zip
python3 scripts/ncbi_datasets_cleanup.py <species>   # -> <species>.chromosomes.fa + <species>_mRNA.gtf
```

`ncbi_datasets_cleanup.py` keeps only full chromosomes, renames them to `chr1..chrN`, and
emits a GTF with `gene_id`/`transcript_id` attributes that the downstream tools expect.

---

## Pipeline steps

Each Snakemake rule runs inside a dedicated container (see [Containers](#containers)).

| # | Rule | Tool | Container | Purpose |
|---|------|------|-----------|---------|
| 1 | `prefetch_sra` *(optional)* | prefetch / fasterq-dump | `sra_tools` | Download RNA-seq runs from SRA and convert to FASTQ |
| 2 | `fastp` | fastp | `fastp` | Adapter/quality trimming of raw reads |
| 3 | `star_genome_generate` | STAR | `star` | Build the STAR genome index from genome + annotation |
| 4 | `star_align` | STAR | `star` | Spliced alignment of trimmed reads → sorted BAM |
| 5 | `samtools_index` | samtools | `samtools` | Index each BAM |
| 6 | `stringtie_assemble` | StringTie | `stringtie` | Reference-guided transcript assembly per sample |
| 7 | `stringtie_merge` | StringTie | `stringtie` | Merge per-sample assemblies into one transcriptome |
| 8 | `extract_candidates` | awk + gffread | `gffread` | Keep novel `MSTRG` transcripts → candidate GTF/FASTA/IDs |
| 9 | `feelnc_filter` | FEELnc_filter.pl | `feelnc` | Remove candidates overlapping mRNA / too short |
| 10 | `cpat` | CPAT (plant model) | `cpat` | Coding-probability score per candidate |
| 11 | `lncfinder` | RNAfold + LncFinder (R) | `lncfinder` | Secondary-structure-aware coding/non-coding SVM |
| 12 | `lncboost` | PlantLncBoost (CatBoost) | `lncboost` | Plant lncRNA classifier |
| 13 | `diamond_makedb` | DIAMOND | `diamond` | Build the protein database (once) |
| 14 | `diamond_blastx` | DIAMOND | `diamond` | Detect protein-coding candidates (translated hits) |
| 15 | `intersect_predictions` | R (base) | `intersect` | Intersect the four tools, subtract protein hits |
| 16 | `final_lncRNA_gtf` | python | `python_utils` | Filter candidate GTF to the consensus IDs |
| 17 | `final_lncRNA_fasta` | gffread | `gffread` | Extract consensus lncRNA sequences |
| 18 | `feelnc_classifier` | FEELnc_classifier.pl | `feelnc` | Classify lncRNAs (intergenic/antisense, etc.) |
| 19 | `truth_cdhit` *(optional)* | CD-HIT | `cd_hit` | Cluster/condense the truth lncRNA set |
| 20 | `minimap2_truth_alignment` *(optional)* | minimap2 | `minimap2` | Benchmark predictions vs the truth set (TP/FP/FN) |

**The consensus logic** (rule 15, `workflow/scripts/intersect_ids.R`): a candidate is a
high-confidence lncRNA only if it is called *non-coding* by **FEELnc, CPAT, LncFinder and
PlantLncBoost simultaneously** and has **no significant DIAMOND protein hit**
(`evalue ≤ 1e-5`). Per-tool cutoffs (`cpat_coding_probability_cutoff`, `lncboost_threshold`,
`diamond_evalue`, …) are set in `config/config.yaml`.

---

## Containers

15 single-tool images, built from the definition files in `containers/*.def`. They are **not**
stored in git (≈3 GB); build them with `scripts/build_singularity_images.sh`. See
[`docs/CONTAINERS.md`](docs/CONTAINERS.md) for the full per-image breakdown.

| Image (`.sif`) | Base | Software (pinned) | Used by |
|----------------|------|-------------------|---------|
| `sra_tools` | micromamba | sra-tools 3.1.1, pigz | prefetch_sra |
| `fastp` | micromamba | fastp 0.23.4 | fastp |
| `star` | micromamba | STAR 2.7.11b | star_genome_generate, star_align |
| `samtools` | micromamba | samtools 1.20 | samtools_index |
| `stringtie` | micromamba | StringTie 2.2.1 | stringtie_assemble/merge |
| `gffread` | micromamba | gffread 0.12.9 | extract_candidates, final_lncRNA_fasta |
| `feelnc` | ubuntu 22.04 | FEELnc (+ BioPerl) | feelnc_filter, feelnc_classifier |
| `cpat_plant` | ubuntu 18.04 | CPAT 1.2.4 + Plant-LncRNA model | cpat |
| `plant_lnc_boost` | python 3.9 | PlantLncBoost (CatBoost) + model | lncboost |
| `lncfinder` | ubuntu 22.04 | R + LncFinder + ViennaRNA 2.7.2 (RNAfold) | lncfinder |
| `diamond` | micromamba | DIAMOND 2.1.8 | diamond_makedb/blastx |
| `intersect_ids` | rocker/r-ver 4.3.3 | base R | intersect_predictions |
| `python_utils` | python 3.11-slim | python 3 | final_lncRNA_gtf |
| `minimap2` | micromamba | minimap2 2.28 | minimap2_truth_alignment |
| `cd_hit` | micromamba | cd-hit 4.8.1 | truth_cdhit |

The `cpat_plant` and `plant_lnc_boost` images clone
[Plant-LncRNA-pipeline-v2](https://github.com/xuechantian/Plant-LncRNA-pipeline-v2) at build
time to bundle the plant hexamer table, logit model and CatBoost model (build node needs
internet).

---

## Outputs

Written to `{species}/03_outputs/10_final/`:

| File | Description |
|------|-------------|
| `{prefix}_final_ids.txt` | Consensus high-confidence lncRNA transcript IDs |
| `{prefix}_lncRNA.gtf` | GTF of the consensus lncRNAs |
| `{prefix}_lncRNA.fasta` | FASTA sequences of the consensus lncRNAs |
| `{prefix}_lncRNA_classes.txt` | FEELnc classification (intergenic / antisense / …) |
| `{prefix}_intersection_summary.csv` | Per-candidate pass/fail table for every tool |

Optional validation (`run_truth_alignment: true`) adds
`{species}/03_outputs/11_verify/{prefix}_truth_vs_pred.metrics.txt` with TP/FP/FN counts.

---

## How to run

```bash
# 1. Environment (Snakemake)
conda env create -f workflow/envs/snakemake.yaml
conda activate iwgc-lnc-snakemake

# 2. Build the containers (build node needs internet; ≈3 GB total)
scripts/build_singularity_images.sh            # or: scripts/build_singularity_images.sh --remote

# 3. Edit config/config.yaml for your species (paths, prefix, samples, cutoffs)

# 4. Dry-run, then run with Singularity
snakemake --snakefile workflow/Snakefile --configfile config/config.yaml --cores 32 --dry-run
snakemake --snakefile workflow/Snakefile --configfile config/config.yaml --cores 32 --use-singularity
```

**On a Slurm cluster**, wrap the `snakemake` call in an `sbatch` script and run it on a
compute node (see `scripts/run_pipeline_ath.sb` for the Arabidopsis example). Note: the
`apptainer`/`singularity` binary must be on `PATH` inside the job, and bind the project tree
with `APPTAINER_BINDPATH`.

A quick **smoke test** validates the custom R/Python logic without containers:

```bash
scripts/run_smoke_test.sh
```

---

## Test results (Arabidopsis TAIR10)

End-to-end validation run on the MSU HPCC. See [`docs/TEST_RESULTS.md`](docs/TEST_RESULTS.md)
for the full write-up (including four bugs found and fixed).

**Inputs used**

| Input | Value |
|-------|-------|
| Genome / annotation | TAIR10 (RefSeq `GCF_000001735.4`) — 5 chromosomes, 55,937 mRNAs |
| RNA-seq | `SRR2073143`, `SRR1688325` — single-end, 50 bp, subsampled to 6M reads each |
| Protein DB | UniProt SwissProt (574,627 sequences) |
| Truth set | none (validation step disabled for this functional test) |

**Pipeline funnel**

```
61,237 merged transcripts → 1,375 MSTRG candidates → 394 scored by ≥1 tool
   FEELnc non-coding   :   9
   CPAT  non-coding    : 286
   LncFinder non-coding: 384
   PlantLncBoost lncRNA: 288
   DIAMOND protein hit : 259
→ 5-way consensus (minus protein hits): 4 high-confidence lncRNAs
```

**Result — 4 high-confidence lncRNAs**

| lncRNA | Chr | Length | FEELnc class |
|--------|-----|--------|--------------|
| `MSTRG.1002.1` | chr1 | 218 bp (2 exons) | — |
| `MSTRG.2278.1` | chr1 | 275 bp (2 exons) | — |
| `MSTRG.9745.1` | chr3 | 585 bp | — |
| `MSTRG.10857.1` | chr3 | 374 bp | antisense / intergenic to `AT3G59260` |

Copies of the actual output files are in
[`docs/test_results/arabidopsis_tair10/`](docs/test_results/arabidopsis_tair10/).

> The small, shallow test set (2 single-end samples) yields only 9 FEELnc survivors, hence a
> small consensus. For real annotation work use many diverse RNA-seq samples at full depth and
> enable the truth-set validation.

---

## Repository layout

```
IWGC_lncRNA_pipeline/
├── README.md                     # this file
├── config/
│   ├── config.yaml               # main configuration
│   └── config.ath.yaml           # Arabidopsis test configuration
├── workflow/
│   ├── Snakefile                 # the pipeline (22 rules)
│   ├── envs/snakemake.yaml       # conda env for Snakemake
│   └── scripts/                  # intersect_ids.R, lncFinder.R, filter_gtf_by_ids.py
├── containers/                   # *.def Singularity definitions (images built into images/)
├── scripts/                      # build script, Slurm templates, cleanup & model files
├── tests/smoke/                  # smoke-test fixtures
├── docs/                         # extended documentation + test results
│   ├── CONTAINERS.md
│   ├── TEST_RESULTS.md
│   └── test_results/arabidopsis_tair10/
└── {species}/                    # per-species inputs (01_inputs) and outputs (03_outputs)
```

---

## Further documentation

- [`docs/CONTAINERS.md`](docs/CONTAINERS.md) — every container image, software, and build notes
- [`docs/TEST_RESULTS.md`](docs/TEST_RESULTS.md) — full Arabidopsis test report + bug fixes
- [`config/config.yaml`](config/config.yaml) — all configurable paths and parameters

# Container images

Every step runs in its own [Singularity/Apptainer](https://apptainer.org/) image, defined by a
`.def` file in [`../containers/`](../containers). One tool per image keeps environments isolated
and the pipeline reproducible.

## Building

```bash
scripts/build_singularity_images.sh                 # all images -> containers/images/
scripts/build_singularity_images.sh hisat2 lncfinder intersect_ids   # a subset
scripts/build_singularity_images.sh --remote        # remote (Sylabs) builder
```

Images (~3.5 GB total) are **git-ignored** — rebuild them from the `.def` files. The build node
needs **internet**: `cpat_plant`, `plant_lnc_boost`, `feelnc`, `lncfinder` and `intersect_ids`
clone upstream repos / download source at build time.

## Image reference

| Image | Base | Key software | Bundled upstream assets | Rules |
|-------|------|--------------|-------------------------|-------|
| `hisat2.sif` | micromamba | HISAT2 2.2.1, samtools 1.20 | — | hisat2_index, hisat2_align |
| `fastp.sif` | micromamba | fastp 0.23.4 | — | fastp |
| `stringtie.sif` | micromamba | StringTie 2.2.1 | — | stringtie_assemble/merge |
| `gffread.sif` | micromamba | gffread 0.12.9 | — | extract_candidates, final_lncRNA_gtf/fasta |
| `diamond.sif` | micromamba | DIAMOND 2.1.8 | — | diamond_makedb/blastx |
| `feelnc.sif` | ubuntu 22.04 | FEELnc (git) + BioPerl | FEELnc scripts | feelnc_filter, feelnc_classifier |
| `cpat_plant.sif` | ubuntu 18.04 | CPAT 1.2.4 (Python 2) + R 3.4 | Plant_Hexamer.tsv | cpat |
| `plant_lnc_boost.sif` | python 3.9 | CatBoost, scikit-learn, BioPython | PlantLncBoost scripts + `.cb` model + `prediction_insersection.sh` | lncboost |
| `lncfinder.sif` | ubuntu 22.04 | R + LncFinder, seqinr, e1071 | **Plant_model.rda + training FASTAs** (v1 repo) | lncfinder |
| `intersect_ids.sif` | rocker/tidyverse 4.3.3 | R + tidyverse + VennDiagram | `prediction_insersection.sh` | intersect_predictions |
| `sra_tools.sif` | micromamba | sra-tools 3.1.1, pigz | — | prefetch_sra |
| `samtools.sif` | micromamba | samtools 1.20 | — | (helper) |
| `minimap2.sif` | micromamba | minimap2 2.28 | — | (optional truth eval) |
| `cd_hit.sif` | micromamba | cd-hit 4.8.1 | — | (optional truth clustering) |
| `python_utils.sif` | python 3.11-slim | Python 3 | — | (helper) |

## Containers built / changed to match upstream

- **`hisat2.sif` (new)** — upstream aligns with HISAT2 (not STAR). Includes samtools so the rule
  pipes `hisat2 … | samtools sort` into one sorted BAM, identical to `hisat2 -S x.sam; samtools
  view -b x.sam | samtools sort`.
- **`lncfinder.sif` (rebuilt)** — clones the predecessor repo
  `Plant-LncRNA-pipline`, which ships `Model/Plant_model.rda` (the plant SVM) and
  `example_data/training_{mRNA,lncRNA}.fasta` (for `make_frequencies`). LncFinder-plant uses
  `SS.features = FALSE`, so no ViennaRNA/RNAfold is needed.
- **`intersect_ids.sif` (rebuilt)** — `rocker/tidyverse` + `VennDiagram`, so the upstream
  `prediction_insersection.sh` runs verbatim (the repo vendors a copy with one cosmetic Venn
  compatibility line — see README → Fidelity).

## HPC storage tip

`.sif` images are large; on clusters with a tight `$HOME` quota keep them on scratch and point the
`containers:` paths there (e.g. `/mnt/scratch/<user>/iwgc_lncRNA/containers/images/*.sif`), adding
that directory to `APPTAINER_BINDPATH`.

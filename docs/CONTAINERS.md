# Container images

Every step of the pipeline runs inside its own [Singularity/Apptainer](https://apptainer.org/)
image, defined by a `.def` file in [`../containers/`](../containers). One tool per image keeps
software environments isolated and the whole pipeline reproducible.

## Building

```bash
# Build all images into containers/images/
scripts/build_singularity_images.sh

# Build only some
scripts/build_singularity_images.sh fastp star stringtie

# On clusters that require a remote (Sylabs) builder
scripts/build_singularity_images.sh --remote
```

The build script picks `apptainer` (or `singularity`) automatically and writes
`containers/images/<name>.sif`. The resulting images (~3 GB total) are **git-ignored** — they
are build artifacts, fully reproducible from the `.def` files.

> The build node needs **internet access**: the `cpat_plant`, `plant_lnc_boost`, `feelnc` and
> `lncfinder` images clone upstream repositories / download source at build time.

## Image reference

| Image | Definition | Base image | Key software (pinned) | Provides | Rules |
|-------|-----------|------------|-----------------------|----------|-------|
| `sra_tools.sif` | `sra_tools.def` | micromamba 1.5.10 | sra-tools 3.1.1, pigz | `prefetch`, `fasterq-dump` | prefetch_sra |
| `fastp.sif` | `fastp.def` | micromamba 1.5.10 | fastp 0.23.4 | `fastp` | fastp |
| `star.sif` | `star.def` | micromamba 1.5.10 | STAR 2.7.11b | `STAR` | star_genome_generate, star_align |
| `samtools.sif` | `samtools.def` | micromamba 1.5.10 | samtools 1.20 | `samtools` | samtools_index |
| `stringtie.sif` | `stringtie.def` | micromamba 1.5.10 | StringTie 2.2.1 | `stringtie` | stringtie_assemble, stringtie_merge |
| `gffread.sif` | `gffread.def` | micromamba 1.5.10 | gffread 0.12.9 | `gffread` | extract_candidates, final_lncRNA_fasta |
| `diamond.sif` | `diamond.def` | micromamba 1.5.10 | DIAMOND 2.1.8 | `diamond` | diamond_makedb, diamond_blastx |
| `minimap2.sif` | `minimap2.def` | micromamba 1.5.10 | minimap2 2.28 | `minimap2` | minimap2_truth_alignment |
| `cd_hit.sif` | `cd_hit.def` | micromamba 1.5.10 | cd-hit 4.8.1 | `cd-hit` | truth_cdhit |
| `python_utils.sif` | `python_utils.def` | python 3.11-slim | Python 3.11 | `python3` | final_lncRNA_gtf |
| `intersect_ids.sif` | `intersect_ids.def` | rocker/r-ver 4.3.3 | base R 4.3.3 | `Rscript` | intersect_predictions |
| `feelnc.sif` | `feelnc.def` | ubuntu 22.04 | FEELnc (git `tderrien/FEELnc`) + BioPerl | `FEELnc_filter.pl`, `FEELnc_classifier.pl` | feelnc_filter, feelnc_classifier |
| `cpat_plant.sif` | `cpat_plant.def` | ubuntu 18.04 | CPAT 1.2.4 (Python 2) + Plant-LncRNA-pipeline-v2 | `cpat.py` + `Plant_Hexamer.tsv` | cpat |
| `plant_lnc_boost.sif` | `plant_lnc_boost.def` | python 3.9-slim | CatBoost, scikit-learn, BioPython + Plant-LncRNA-pipeline-v2 | `Feature_extraction.py`, `PlantLncBoost_prediction.py` + `.cb` model | lncboost |
| `lncfinder.sif` | `lncfinder.def` | ubuntu 22.04 | R + LncFinder, seqinr, e1071 + ViennaRNA 2.7.2 | `Rscript`, `RNAfold` | lncfinder |

## Notes on the model containers

- **`cpat_plant`** — CPAT 1.2.4 runs under Python 2 (legacy). The image clones
  `Plant-LncRNA-pipeline-v2` into `/opt`, providing the plant-trained hexamer table
  (`/opt/Plant-LncRNA-pipeline-v2/Model/Plant_Hexamer.tsv`). The logit model
  (`scripts/Plant.logit.v2.RData`) is supplied from this repository via `config.cpat_logit`.
- **`plant_lnc_boost`** — bundles the CatBoost model
  (`/opt/Plant-LncRNA-pipeline-v2/Model/PlantLncBoost_model.cb`) and the feature-extraction /
  prediction scripts under `/opt/Plant-LncRNA-pipeline-v2/PlantLncBoost/Script/`.
- **`lncfinder`** — compiles ViennaRNA 2.7.2 from source to provide `RNAfold`, plus the R
  packages `LncFinder`, `seqinr`, `e1071`. The pipeline uses the built-in `wheat` plant SVM
  model with secondary-structure features.
- **`feelnc`** — installs the BioPerl stack required by FEELnc and symlinks
  `FEELnc_filter.pl` / `FEELnc_classifier.pl` onto `PATH`.

## Storage tip (HPC)

The `.sif` images are large. On clusters with a tight `$HOME` quota, keep them on scratch and
point the `containers:` paths in the config there, e.g.
`/mnt/scratch/<user>/iwgc_lncRNA/containers/images/*.sif`, and add that directory to
`APPTAINER_BINDPATH`.

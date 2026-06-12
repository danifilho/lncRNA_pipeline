# Test results — Arabidopsis thaliana (TAIR10)

End-to-end validation of the faithful, containerized pipeline on the **MSU HPCC (ICER)**
(Snakemake + Apptainer). Confirms the workflow runs to completion and produces a non-empty
high-confidence lncRNA set, command-for-command equivalent to Plant-LncRNA-pipeline-v2.

## Inputs

| Input | Value | Source |
|-------|-------|--------|
| Genome + annotation | TAIR10, RefSeq `GCF_000001735.4` (5 chromosomes, 55,937 mRNAs) | NCBI |
| RNA-seq | `SRR2073143`, `SRR1688325` — single-end, 50 bp, 6,000,000 reads each | ENA |
| Protein DB | UniProt SwissProt (574,627 sequences) | UniProt |

Config: [`config/config.ath.yaml`](../config/config.ath.yaml) (species `A.thaliana`, prefix
`ath`). Slurm wrapper: [`scripts/run_pipeline_ath.sb`](../scripts/run_pipeline_ath.sb).

## Funnel

| Stage | Count |
|-------|------:|
| `MSTRG` candidate transcripts | 2,674 |
| FEELnc candidate lncRNAs (`candidate_lncRNA.txt`) | 14 |
| LncFinder-plant non-coding | 585 |
| PlantLncBoost lncRNA (label = 1) | 593 |
| CPAT non-coding (real `coding_prob<0.46`) | 582 |
| **Consensus (`prediction_insersection.sh`, verbatim)** | **10** |

## Result — 10 high-confidence lncRNAs

```
MSTRG.1000.1  MSTRG.2055.1  MSTRG.2326.2  MSTRG.3725.1  MSTRG.5145.1
MSTRG.9727.1  MSTRG.10838.1 MSTRG.11740.1 MSTRG.14179.1 MSTRG.16769.1
```

Output files in [`test_results/arabidopsis_tair10/`](test_results/arabidopsis_tair10/):
`Final_lncRNA_results.txt`, `lncRNA.gtf`, `lncRNA.fasta`, `lncRNA_classes.txt`,
`Venn_pred_lncRNA.pdf`.

## Reproducing

```bash
sbatch scripts/run_pipeline_ath.sb   # one Slurm job: Snakemake + Apptainer, ~20 min
```

## Notes / caveats

- **CPAT is effectively permissive** in the upstream `prediction_insersection.sh` (see
  [README → Fidelity](../README.md#fidelity-to-upstream)). Reproduced verbatim; the consensus is
  driven by FEELnc ∩ PlantLncBoost ∩ LncFinder ∩ ¬protein.
- The test set is intentionally small (2 single-end samples) — only 14 FEELnc survivors, hence a
  small consensus. This validates correctness, not biological completeness. For real annotation
  use many diverse RNA-seq samples at full depth.
- Outputs for the HPCC test run are written under `/mnt/scratch/dasilvaf/iwgc_lncRNA/` to avoid
  the tight `$HOME` quota; container images live on scratch as well.

## How this matches upstream

Aligner (HISAT2), assembler (StringTie), all four classifiers (FEELnc, CPAT-plant with
`Plant.logit`, LncFinder-plant with `Plant_model.rda` + `SS.features=FALSE`, PlantLncBoost),
DIAMOND, and the intersection script are the upstream tools/commands/models. The three packaging
deviations (logit serialization format, training-data source, Venn compat line) are documented in
the README and change neither the method nor the result.

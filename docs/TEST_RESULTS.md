# Test results — Arabidopsis thaliana (TAIR10)

End-to-end validation of the containerized pipeline on the **MSU HPCC (ICER)**
(Snakemake). Confirms the workflow runs to completion and produces a non-empty
lncRNA set, command-for-command equivalent to Plant-LncRNA-pipeline-v2.

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

---

# Test results — Glycine max (Phytozome Wm82.a6.v1)

The upstream README's own example data (Plant-LncRNA-pipeline-v2 §4), run end-to-end on the MSU
HPCC. Slurm job `9944833`, `COMPLETED`, 1 h 54 m on 24 cores / 96 GB.

## Inputs

| Input | Value | Source |
|-------|-------|--------|
| Genome + annotation | **Glycine max Wm82.a6.v1** — `Gmax_880_v6.0.fa` + `Gmax_880_Wm82.a6.v1.gene_exons.gff3` (→ GTF via `gffread -T`) | Phytozome 13 (JGI) |
| RNA-seq | `SRR1174214`, `SRR1174217`, `SRR1174218`, `SRR1174232` — **paired-end**, ~25–37 M pairs each | SRA / ENA |
| Library | strand-specific (`strand_specific: true` → `hisat2 --rna-strandness RF`, `stringtie --rf`) | as in the upstream README |
| Protein DB | UniProt SwissProt | UniProt |

Config: [`config/config.gmax.yaml`](../config/config.gmax.yaml). Slurm wrapper:
[`scripts/run_pipeline_gmax.sb`](../scripts/run_pipeline_gmax.sb). Step-by-step:
[`docs/RUNNING_GMAX_EXAMPLE.md`](RUNNING_GMAX_EXAMPLE.md).

## Funnel

| Stage | Count |
|-------|------:|
| `MSTRG` candidate transcripts | 75,327 |
| FEELnc candidate lncRNAs | 53,515 |
| LncFinder-plant non-coding | 62,564 |
| PlantLncBoost lncRNA (label = 1) | 62,194 |
| DIAMOND protein hits removed (`pident>60 & e<1e-5`, unique) | 36,164 |
| **Consensus (`prediction_insersection.sh`, verbatim)** | **23,606** |

## Classification of the 23,606 lncRNAs (FEELnc_classifier, upstream §9)

| Category | Count |
|----------|------:|
| antisense / exonic | 21,341 |
| intergenic | 1,153 |
| downstream | 231 |
| upstream | 199 |
| bidirectional | 191 |
| intronic | 101 |

Light outputs archived in [`test_results/gmax_wm82/`](test_results/gmax_wm82/)
(`Final_lncRNA_results.txt`, the category files, `Venn_pred_lncRNA.pdf`). The full
`lncRNA.gtf` (14 MB) / `lncRNA.fasta` (27 MB) / `lncRNA_classes.txt` (15 MB) stay on scratch.

## Notes

- **Why so many vs Arabidopsis's 10?**

> **Is 23,606 reasonable?**  See
> [RESULTS_INTERPRETATION.md](RESULTS_INTERPRETATION.md) for the full comparison and the analysis of
> the upstream CPAT parsing issue (CPAT effectively does not filter).

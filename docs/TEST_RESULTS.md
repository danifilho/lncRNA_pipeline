# Test results — Arabidopsis thaliana (TAIR10)

End-to-end validation of the pipeline on the **MSU HPCC (ICER)** using Snakemake +
Apptainer. This run confirmed the workflow runs to completion and produces a non-empty
high-confidence lncRNA set, and surfaced four real bugs (all fixed — see below).

## Inputs

| Input | Value | Source |
|-------|-------|--------|
| Genome + annotation | TAIR10, RefSeq `GCF_000001735.4` — 5 chromosomes, 55,937 mRNA transcripts | NCBI |
| RNA-seq | `SRR2073143`, `SRR1688325` — single-end, 50 bp, subsampled to 6,000,000 reads each | ENA |
| Protein database | UniProt SwissProt — 574,627 sequences | UniProt |
| Truth set | none (`run_truth_alignment: false`) | — |

Configuration: [`config/config.ath.yaml`](../config/config.ath.yaml)
(species `A.thaliana`, prefix `ath`). Slurm wrapper:
[`scripts/run_pipeline_ath.sb`](../scripts/run_pipeline_ath.sb).

## Pipeline funnel

| Stage | Count |
|-------|------:|
| Merged StringTie transcripts | 61,237 |
| `MSTRG` candidate transcripts (gffread) | 1,375 |
| Candidates scored by ≥1 tool (summary universe) | 394 |
| → FEELnc non-coding | 9 |
| → CPAT non-coding (`coding_prob < 0.46`) | 286 |
| → LncFinder non-coding | 384 |
| → PlantLncBoost lncRNA | 288 |
| → DIAMOND protein hit (`e ≤ 1e-5`) | 259 |
| **5-way consensus (non-coding by all 4, no protein hit)** | **4** |

## Result — 4 high-confidence lncRNAs

| lncRNA | Chr | Strand | Length | Exons | FEELnc class |
|--------|-----|--------|--------|-------|--------------|
| `MSTRG.1002.1`  | chr1 | − | 218 bp | 2 | — |
| `MSTRG.2278.1`  | chr1 | + | 275 bp | 2 | — |
| `MSTRG.9745.1`  | chr3 | − | 585 bp | 1 | — |
| `MSTRG.10857.1` | chr3 | − | 374 bp | 1 | antisense / intergenic to `AT3G59260` (298 bp, convergent) |

Output files (copied into [`test_results/arabidopsis_tair10/`](test_results/arabidopsis_tair10/)):
`ath_final_ids.txt`, `ath_lncRNA.gtf`, `ath_lncRNA.fasta`, `ath_lncRNA_classes.txt`,
`ath_intersection_summary.csv`.

## Bugs found and fixed

The pipeline did not work out of the box — the LncFinder step crashed and, once it ran, the
consensus was silently always empty. Four bugs were fixed:

1. **STAR GTF tag mismatch** (`workflow/Snakefile`, `star_genome_generate`)
   Used `--sjdbGTFtagExonParentTranscript Parent` (a GFF3 convention), but
   `ncbi_datasets_cleanup.py` emits a standard GTF with `transcript_id`. Changed to
   `transcript_id`.

2. **`lncFinder.R` — wrong LncFinder API + case sensitivity**
   - The installed `LncFinder::lnc_finder()` has no `Structures=` argument; it expects a
     `run_RNAfold()` / `read_SS()`-style data frame (`format = "SS"`). Rewrote the script to
     build that structure from the RNAfold output.
   - Candidate transcripts with NA features (no ORF) silently dropped during `predict.svm()`
     and broke the internal `cbind`. Now pre-filtered with `complete.cases(extract_features())`.
   - The gffread candidate FASTA is soft-masked (mixed case) and LncFinder's k-mer features are
     case-sensitive; without `tolower()` every prediction collapsed to ~0.57 "Coding". Fixed.

3. **`intersect_ids.R` — CPAT ID parsing** (`get_cpat_ids`)
   `cpat.py` writes the sequence ID as an *unnamed* first column, so `read.table(header=TRUE)`
   put the IDs into row names and the code read `mRNA_size` numbers as IDs. CPAT IDs therefore
   never matched the `MSTRG.x.y` IDs from the other tools, so **no transcript ever passed
   CPAT** and the consensus was always empty. Fixed to use row names.

4. **`intersect_ids.R` — LncFinder excluded from the consensus** (`get_lncfinder_ids`)
   The LncFinder results are a *quoted* CSV, but the reader used `quote = ""`, leaving literal
   quotes on the IDs/labels so nothing matched; combined with the rule that drops empty tool
   sets, LncFinder was silently left out of the intersection. Fixed to parse with `read.csv()`.

The smoke test (`scripts/run_smoke_test.sh`) passes after all changes. Bugs (3) and (4) were
the actual reason the final lncRNA set was empty; fixing them produced the 4 lncRNAs above.

## Reproducing

```bash
# inputs prepared under A.thaliana/01_inputs and resources/uniprot (see config.ath.yaml)
sbatch scripts/run_pipeline_ath.sb        # single Slurm job, Snakemake + Apptainer
```

## Caveats / next steps

- The test set is intentionally tiny (2 single-end samples, 6M reads each) → only 9 FEELnc
  survivors → small consensus. This validates correctness, not biological completeness.
- For real annotation work: use many diverse RNA-seq samples at full depth (ideally the full
  PLncDB set for the species), and enable `run_truth_alignment: true` with a truth lncRNA set
  to measure sensitivity/precision.

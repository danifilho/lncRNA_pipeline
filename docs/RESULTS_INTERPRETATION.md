# Interpreting the results & next steps

This note explains how to read the consensus counts produced by the pipeline (in particular the
Glycine max run), documents a parsing issue inherited from the upstream intersection script, and
lists concrete next steps. Nothing here has been applied to the faithful pipeline — these are
recommendations.

## How many lncRNAs should we expect? (Glycine max)

| Reference | lncRNAs | Notes |
|-----------|--------:|-------|
| Typical single genome-wide soybean study | ~1,500–1,600 | few samples, stringent (e.g. 1,563; 1,578 with 1,010 intergenic / 465 antisense / 84 intronic) |
| **PLncDB v2** (curated DB, many studies combined) | **~11,124 unique** (12,577 records) | the database this pipeline's truth step compares against |
| **This pipeline** (4 deep runs, Wm82.a6.v1) | **23,606** | faithful run, job 9944833 |

The 23,606 is the **right order of magnitude** (thousands–tens of thousands) but **~2× the entire
PLncDB** and **~15× a typical single study** → it is an **over-estimate**. This is **not a pipeline
failure** — it runs the upstream method faithfully — it reflects the test setup plus known looseness
of the upstream approach.

## Why the count is inflated

1. **CPAT effectively does not filter** (see next section) — the consensus is really
   FEELnc ∩ LncFinder ∩ PlantLncBoost ∩ ¬protein instead of a true 4-tool intersection.
2. **Isoforms are not collapsed** — counts are transcript-level (`MSTRG.x.y`), not gene-level.
3. **FEELnc kept 71 %** of candidates (53,515 / 75,327) — this depends on the quality of the
   `gffread`-converted Phytozome annotation; a longest-isoform reference would be stricter.
4. **Classification is ~90 % antisense/exonic** (21,341 of 23,606), whereas real studies are
   dominated by *intergenic* lncRNAs — a sign of over-calling (likely spurious antisense
   transcription captured by the RF-stranded libraries).

## The CPAT parsing issue (upstream intersection script)

`workflow/scripts/prediction_insersection.sh` is vendored **verbatim** from
Plant-LncRNA-pipeline-v2; its CPAT line is:

```r
CPAT <- read_delim(cpat_output, delim = "\t") %>% filter(coding_prob < 0.46) %>% pull(mRNA_size)
```

**CPAT 1.2.4 writes the sequence ID as an unnamed first column** — the header has **5** names but
each data row has **6** fields:

```
header:  mRNA_size   ORF_size   Fickett_score   Hexamer_score   coding_prob          (5 names)
data  :  <ID>        <mRNA_size> <ORF_size>     <Fickett>       <Hexamer> <coding_prob>  (6 fields)
```

`read_delim` assigns the 5 names to the **first 5 columns**, so:

- the column named `mRNA_size` is actually the **ID** → `pull(mRNA_size)` returns the right IDs ✓
- the column named `coding_prob` is **not** the coding probability — the real `coding_prob`
  (the 6th field) overflows past the 5 names. So **`filter(coding_prob < 0.46)` never filters on the
  coding probability.**

Verified with **readr 2.1.5 / tidyverse 2.0.0** (the version in `intersect_ids.sif`) on the soybean
CPAT output:

| Expression | IDs kept |
|------------|---------:|
| upstream `filter(coding_prob<0.46) %>% pull(mRNA_size)` | **75,291 / 75,327** (~all) |
| correct parse (real `coding_prob < 0.46`) | **61,635** |

So CPAT contributes essentially nothing to the consensus.

### Is this present in the original repo? — what we can and cannot claim

- **Certain:** the script is verbatim upstream, and the column misalignment is a **deterministic**
  consequence of CPAT 1.2.4's format (6 data fields, 5 header names). The real `coding_prob` is the
  6th field and is **never** named `coding_prob`, in **any** readr version — so the filter does not
  use the coding probability regardless of versions.
- **Version-dependent (uncertain):** the **exact** number kept. readr ≥ 2 merges the overflow into
  the last column (a string), so `coding_prob < 0.46` becomes a string comparison that matches
  almost everything; readr 1.x would instead leave `coding_prob` holding the **Hexamer** column
  (also wrong, but a different result). The upstream repo **does not pin package versions**, so we
  cannot know which behavior the original authors obtained, nor whether they noticed it.

We therefore reproduce the **published script's behavior with current tidyverse, by choice**
(keeping it verbatim). We do **not** claim the authors' own published results are wrong.

## Validation against PLncDB (minimap2)

We aligned the 23,606 predicted lncRNAs against known soybean lncRNAs from **PLncDB** (5,767
non-redundant *Glycine max* sequences, taxid 3847, via RNAcentral) with `minimap2 -x cdna`
(`scripts/validate_vs_truth.sh`; full metrics in
[`test_results/gmax_wm82/validation_vs_PLncDB.txt`](test_results/gmax_wm82/validation_vs_PLncDB.txt)).

| identity cutoff | predictions matching a known lncRNA | known lncRNAs recovered (recall) |
|---|---:|---:|
| 0.80 | 764 (3.2 %) | 619 / 5,767 (**10.7 %**) |
| 0.90 | 511 (2.2 %) | 455 / 5,767 (**7.9 %**) |
| 0.95 | 394 (1.7 %) | 373 / 5,767 (**6.5 %**) |

**Reading these numbers honestly:**

- **Recall is low (~7–11 %).** We recover only a small fraction of known soybean lncRNAs. This is
  expected: PLncDB aggregates many studies and tissues, whereas this run used only **4 RNA-seq
  libraries (mostly root)**, and lncRNAs are highly tissue/condition-specific — a small sample set
  captures few of them.
- **Only ~2–3 % of predictions match a known lncRNA.** A low match rate alone does **not** prove the
  rest are false positives (PLncDB is not exhaustive; genuinely novel lncRNAs are possible). **But
  combined with the inflation signals** (CPAT not filtering, ~90 % antisense, no isoform collapse),
  it indicates that a large fraction of the 23,606 are low-confidence / spurious rather than
  bona-fide lncRNAs.
- Assemblies differ (predictions Wm82.a6 / Phytozome vs PLncDB Wm82.a4 / NCBI), but the same cultivar
  is ~identical between versions, so this low overlap is **real, not assembly noise**.

**Takeaway:** treat the 23,606 as a *permissive candidate set*, not a curated annotation. The low
recall is driven mostly by limited input data; the low match rate plus the inflation factors mean the
next steps below would materially tighten the result.

## Next steps / possible solutions (not applied — would deviate from "verbatim upstream")

1. **Fix the CPAT parse (1 line):** read the CPAT output by row names / explicit column index so
   `coding_prob` is the real value (as done in this repo's own `intersect_ids.R` helper). Expected
   to pull the consensus toward the literature range.
2. **Collapse to one transcript per gene** (longest isoform) for gene-level counts.
3. **Use a longest-isoform reference GTF** for FEELnc to tighten overlap filtering.
4. **Truth-set validation — done** (see section above): minimap2 vs PLncDB gives ~7–11 % recall and
   ~2–3 % of predictions matching known lncRNAs. To raise recall, add many more RNA-seq libraries
   across tissues/conditions; to raise the match rate, apply steps 1–3.
5. **Confirm library strandness** (e.g. RSeQC `infer_experiment.py`) before trusting the
   antisense-heavy classification.

## References

- PLncDB v2 — plant lncRNA database (truth set; ~11k soybean lncRNAs)
- Tian et al., *PlantLncBoost*, New Phytologist (2025), doi:10.1111/nph.70211
- Tian et al., *Plant-LncPipe*, Horticulture Research 11(4):uhae041 (2024)
- Golicz et al., soybean lncRNAs & small peptides, Plant Physiol. 182(3):1359 (single-study ~1.5k)

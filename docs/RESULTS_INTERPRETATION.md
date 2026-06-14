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

## Next steps / possible solutions (not applied — would deviate from "verbatim upstream")

1. **Fix the CPAT parse (1 line):** read the CPAT output by row names / explicit column index so
   `coding_prob` is the real value (as done in this repo's own `intersect_ids.R` helper). Expected
   to pull the consensus toward the literature range.
2. **Collapse to one transcript per gene** (longest isoform) for gene-level counts.
3. **Use a longest-isoform reference GTF** for FEELnc to tighten overlap filtering.
4. **Enable the truth-set validation** (`run_truth_alignment`: cd-hit + minimap2 vs the PLncDB
   soybean lncRNAs) — already in the pipeline, currently off — to report TP/FP/FN against the ~11k
   PLncDB set and turn the count into a precision/recall statement.
5. **Confirm library strandness** (e.g. RSeQC `infer_experiment.py`) before trusting the
   antisense-heavy classification.

## References

- PLncDB v2 — plant lncRNA database (truth set; ~11k soybean lncRNAs)
- Tian et al., *PlantLncBoost*, New Phytologist (2025), doi:10.1111/nph.70211
- Tian et al., *Plant-LncPipe*, Horticulture Research 11(4):uhae041 (2024)
- Golicz et al., soybean lncRNAs & small peptides, Plant Physiol. 182(3):1359 (single-study ~1.5k)

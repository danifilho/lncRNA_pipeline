# Strict-variant experiment (Glycine max)

**Self-contained experiment — does NOT modify the faithful pipeline.** It post-processes a
*finished* g.max run (no re-alignment/assembly) to test a hypothesis, and keeps all originals
(`workflow/`, `config/`, `docs/`, the verbatim intersection) untouched.

## Hypothesis

The verbatim run yields **23,606** lncRNAs — an over-estimate (≈2× PLncDB, ≈15× a single study;
see [`../docs/RESULTS_INTERPRETATION.md`](../docs/RESULTS_INTERPRETATION.md)). Two changes were
proposed as likely fixes:

1. **Fix the CPAT parse** so CPAT filters on the *real* `coding_prob` (< 0.46), instead of the
   upstream `read_delim` column-misalignment that makes CPAT match ~everything.
2. **Collapse isoforms** to one (longest) transcript per StringTie gene (`MSTRG.x`).

**Question:** do these tighten the result and improve agreement with PLncDB?

## What it does / how to run

- [`strict_consensus.py`](strict_consensus.py) — reads the finished run's per-tool outputs, parses
  CPAT correctly, recomputes `FEELnc ∩ CPAT ∩ LncFinder ∩ PlantLncBoost − protein`, and collapses to
  the longest isoform per gene.
- [`run_strict_variant.sh`](run_strict_variant.sh) — runs the above, extracts the gene-representative
  sequences, and re-validates against PLncDB with `../scripts/validate_vs_truth.sh` (minimap2).

```bash
bash strict_variant/run_strict_variant.sh   # reads scratch outputs; writes strict_variant/results/
```

## Results

Per-tool positive sets (unchanged candidate pool of 75,327 transcripts):

| Tool | positives |
|------|----------:|
| FEELnc candidates | 53,515 |
| CPAT non-coding (**real** `coding_prob < 0.46`) | 61,635 |
| LncFinder non-coding | 62,564 |
| PlantLncBoost lncRNA | 62,194 |
| DIAMOND protein hits removed | 36,164 |

| Consensus | count |
|-----------|------:|
| Verbatim (CPAT not filtering) — transcripts | **23,606** |
| Strict (real CPAT filter) — transcripts | **23,411**  (−195, **−0.8 %**) |
| Strict + isoform collapse — genes (longest isoform) | **18,575**  (−21 % vs strict transcripts) |

Validation vs PLncDB (5,767 known soybean lncRNAs, `minimap2 -x cdna`):

| Set | recall @0.80 / 0.90 / 0.95 | predictions matching known @0.80 / 0.90 / 0.95 |
|-----|----------------------------|-----------------------------------------------|
| Verbatim (23,606 transcripts) | 10.7 % / 7.9 % / 6.5 % | 3.2 % / 2.2 % / 1.7 % |
| **Strict (18,575 genes)** | **10.2 % / 7.6 % / 6.2 %** | **3.1 % / 2.1 % / 1.6 %** |

## Conclusion — hypothesis **not supported**

- **Fixing the CPAT parse removes only 195 transcripts (−0.8 %).** The bug is real, but it has almost
  no effect on *this* consensus: FEELnc ∩ LncFinder ∩ PlantLncBoost already exclude nearly all
  coding transcripts, so a correctly-filtering CPAT agrees with them on all but 195. (Our earlier
  write-up over-emphasised CPAT as the inflation driver — this experiment corrects that.)
- **Isoform collapse gives a modest −21 %** (23,411 → 18,575 genes) but the number is still far
  above PLncDB (~5.8 k RNAcentral / ~11 k full) and a single study (~1.5 k).
- **Validation is essentially unchanged** — recall and match-rate move by <1 point. The strict set is
  not meaningfully "closer" to known lncRNAs.

**So the over-estimation is driven by the candidate pool itself** — 75 k assembled transcripts from
4 libraries, with the three plant tools agreeing on ~23 k as non-coding even though only ~2–3 %
correspond to known lncRNAs — **not by the CPAT quirk or isoforms.** Tightening it would require
deeper/more-diverse RNA-seq (more tissues → better assembly + recall), an **expression filter**
(e.g. a TPM cutoff), a **stricter candidate definition** (drop single-exon / low-coverage
transcripts), and/or higher per-tool thresholds — not the two tweaks tested here.

## Files

- `results/strict_final_transcript_ids.txt` — strict consensus transcripts (23,411)
- `results/strict_final_gene_ids.txt` — longest isoform per gene (18,575)
- `results/per_tool_counts.txt`, `results/validation_vs_PLncDB.txt`

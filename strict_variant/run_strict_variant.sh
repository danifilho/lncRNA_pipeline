#!/usr/bin/env bash
# Strict-variant experiment: post-process a finished g.max run (no re-alignment).
#   1) strict consensus (real CPAT filter) 2) collapse isoforms 3) revalidate vs PLncDB.
set -euo pipefail
SC=/mnt/scratch/dasilvaf/iwgc_lncRNA
O=$SC/g.max/03_outputs
REPO=/mnt/ffs24/home/dasilvaf/projects/IWGC_lncRNA_pipeline
OUT=$REPO/strict_variant/results
WORK=$SC/strict_variant
mkdir -p "$OUT" "$WORK"
export PATH="$PATH:/mnt/ffs24/home/dasilvaf/miniforge3/bin"

echo ">> 1. strict consensus + isoform collapse"
python3 "$REPO/strict_variant/strict_consensus.py" \
  --feelnc    "$O/05_FEELnc/candidate_lncRNA.txt" \
  --cpat      "$O/06_cpat/CPAT_plant.output" \
  --lncfinder "$O/08_lncFinder/plant-lncFinder.txt" \
  --boost     "$O/07_lncBoost/PlantLncBoost_prediction.csv" \
  --diamond   "$O/09_diamond/uniprotoutput.txt" \
  --fasta     "$O/04_gffread/candidate_transcript.fasta" \
  --outdir    "$OUT"

echo ">> 2. extract longest-isoform-per-gene sequences"
awk 'NR==FNR{k[$1]=1;next} /^>/{id=substr($1,2);p=(id in k)} p' \
  "$OUT/strict_final_gene_ids.txt" "$O/04_gffread/candidate_transcript.fasta" > "$WORK/strict_lncRNA.fasta"
echo "   strict gene-repr seqs: $(grep -c '^>' "$WORK/strict_lncRNA.fasta")"

echo ">> 3. validate strict set vs PLncDB"
bash "$REPO/scripts/validate_vs_truth.sh" "$WORK/strict_lncRNA.fasta" "$SC/truth/gmax_PLncDB_truth.fa" "$WORK/validation"
cp "$WORK/validation/metrics.txt" "$OUT/validation_vs_PLncDB.txt"
echo ">> done. results in $OUT"

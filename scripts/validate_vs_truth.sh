#!/usr/bin/env bash
# Benchmark predicted lncRNAs against a known truth set (e.g. PLncDB) with minimap2.
# BENCHMARKING ADD-ON -- not part of upstream Plant-LncRNA-pipeline-v2.
#
#   bash scripts/validate_vs_truth.sh <pred.fasta> <truth.fasta> <out_dir>
#
# Truth set used here: PLncDB soybean lncRNAs (taxid 3847) from RNAcentral, already
# non-redundant. If your truth comes from multiple sources, pre-cluster it first with
# cd-hit-EST (nucleotide!), e.g.:
#   cd-hit-est -i truth.fa -o truth_nr.fa -c 0.90 -n 9 -M 8000 -T 8 -l 200
set -euo pipefail

PRED="${1:?pred.fasta}"; TRUTH="${2:?truth.fasta}"; OUT="${3:?out_dir}"
SIF=/mnt/scratch/dasilvaf/iwgc_lncRNA/containers/images
BIND=/mnt/scratch/dasilvaf/iwgc_lncRNA
export PATH="$PATH:/mnt/ffs24/home/dasilvaf/miniforge3/bin"
mkdir -p "$OUT"

echo ">> minimap2 (-x cdna): predictions vs truth"
apptainer exec -B "$BIND" "$SIF/minimap2.sif" \
  minimap2 -x cdna -t 8 "$TRUTH" "$PRED" -o "$OUT/pred_vs_truth.paf" 2>/dev/null

TOTAL_PRED=$(grep -c '^>' "$PRED")
TOTAL_TRUTH=$(grep -c '^>' "$TRUTH")
{
  echo "total_predicted=${TOTAL_PRED}"
  echo "total_truth=${TOTAL_TRUTH}"
  echo "# per-block identity = PAF matches/blocklen (col10/col11)"
  for c in 0.80 0.90 0.95; do
    MP=$(awk -v c="$c" '$11>0 && $10/$11>c{print $1}' "$OUT/pred_vs_truth.paf" | sort -u | wc -l)
    HT=$(awk -v c="$c" '$11>0 && $10/$11>c{print $6}' "$OUT/pred_vs_truth.paf" | sort -u | wc -l)
    awk -v c="$c" -v mp="$MP" -v ht="$HT" -v tp="$TOTAL_PRED" -v tt="$TOTAL_TRUTH" 'BEGIN{
      printf "cutoff=%s  pred_matching_known=%d (%.1f%% of predictions)  known_recovered=%d (%.1f%% recall)\n",
             c, mp, 100*mp/tp, ht, 100*ht/tt}'
  done
} | tee "$OUT/metrics.txt"

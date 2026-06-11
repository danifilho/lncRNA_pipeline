#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT_DIR="${ROOT_DIR}/tests/smoke/out"

rm -rf "${OUT_DIR}"
mkdir -p "${OUT_DIR}"

Rscript "${ROOT_DIR}/workflow/scripts/intersect_ids.R" \
  --feelnc "${ROOT_DIR}/tests/smoke/feelnc.gtf" \
  --cpat "${ROOT_DIR}/tests/smoke/cpat.output" \
  --boost "${ROOT_DIR}/tests/smoke/lncboost_predictions.txt" \
  --lncfinder "${ROOT_DIR}/tests/smoke/lncfinder_results.csv" \
  --diamond "${ROOT_DIR}/tests/smoke/diamond.tsv" \
  --out "${OUT_DIR}/smoke"

diff -u "${ROOT_DIR}/tests/smoke/expected_final_ids.txt" "${OUT_DIR}/smoke_final_ids.txt"

python3 -m py_compile "${ROOT_DIR}/workflow/scripts/filter_gtf_by_ids.py"
python3 "${ROOT_DIR}/workflow/scripts/filter_gtf_by_ids.py" \
  --ids "${OUT_DIR}/smoke_final_ids.txt" \
  --gtf "${ROOT_DIR}/tests/smoke/feelnc.gtf" \
  --output "${OUT_DIR}/smoke_filtered.gtf"

grep -q 'transcript_id "MSTRG.1.1"' "${OUT_DIR}/smoke_filtered.gtf"
if grep -q 'transcript_id "MSTRG.2.1"' "${OUT_DIR}/smoke_filtered.gtf"; then
  echo "Filtered GTF retained a transcript that should have been removed." >&2
  exit 1
fi

if [ -x "${ROOT_DIR}/.snakemake-test/bin/snakemake" ]; then
  SNAKEMAKE="${ROOT_DIR}/.snakemake-test/bin/snakemake"
elif command -v snakemake >/dev/null 2>&1; then
  SNAKEMAKE="$(command -v snakemake)"
else
  SNAKEMAKE=""
fi

if [ -n "${SNAKEMAKE}" ]; then
  mkdir -p "${ROOT_DIR}/.home" "${ROOT_DIR}/.cache"
  HOME="${ROOT_DIR}/.home" XDG_CACHE_HOME="${ROOT_DIR}/.cache" "${SNAKEMAKE}" \
    --snakefile "${ROOT_DIR}/workflow/Snakefile" \
    --configfile "${ROOT_DIR}/tests/smoke/config.yaml" \
    --cores 1 \
    --dry-run \
    --quiet
else
  echo "snakemake is not installed; skipped Snakefile dry-run."
fi

echo "Smoke test passed."

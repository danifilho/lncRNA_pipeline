#!/usr/bin/env bash
# Download paired/single FASTQ from ENA for every accession in a list.
# RUN ON THE DEV/LOGIN NODE (compute nodes have no internet).
#   bash scripts/download_sra_fastq.sh <sra_list.txt> <out_dir>
set -euo pipefail

LIST="${1:?usage: download_sra_fastq.sh <sra_list.txt> <out_dir>}"
OUT="${2:?usage: download_sra_fastq.sh <sra_list.txt> <out_dir>}"
mkdir -p "$OUT"

while read -r acc; do
  [ -z "$acc" ] && continue
  case "$acc" in \#*) continue;; esac
  echo ">> $acc"
  # ask ENA for the fastq ftp paths
  urls=$(curl -s --max-time 60 \
    "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${acc}&result=read_run&fields=fastq_ftp&format=tsv" \
    | tail -1 | cut -f1)
  if [ -z "$urls" ]; then echo "   !! no ENA fastq for $acc (skipping)"; continue; fi
  IFS=';' read -ra parts <<< "$urls"
  for u in "${parts[@]}"; do
    f="$OUT/$(basename "$u")"
    if [ -s "$f" ]; then echo "   have $(basename "$u")"; continue; fi
    echo "   wget $(basename "$u")"
    wget -q -O "$f" "https://${u}"
  done
done < "$LIST"
echo "All downloads finished -> $OUT"

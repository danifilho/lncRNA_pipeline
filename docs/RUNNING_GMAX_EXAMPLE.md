# Running the Glycine max example (Plant-LncRNA-pipeline-v2 §4) on the MSU HPCC

Step-by-step guide to run the faithful pipeline on the upstream README's example data
(Glycine max, Phytozome **Wm82.a6.v1** + the `SRR1174*` RNA-seq runs) and follow it yourself.

Everything is already prepared except the two inputs that need a manual download:
the **Phytozome genome** (JGI login) and the **RNA-seq FASTQ** (dev node, no internet on compute
nodes). Ready-made files: `config/config.gmax.yaml`, `scripts/run_pipeline_gmax.sb`,
`scripts/download_sra_fastq.sh`.

```
ssh amd24
cd /mnt/home/dasilvaf/projects/IWGC_lncRNA_pipeline
```

> Heads-up on scale: the `SRR1174*` runs are **paired-end, ~5 GB each**, aligned to a ~1 Gb
> soybean genome. Expect a multi-hour Slurm job and tens of GB on scratch. Start with a few
> accessions. All large I/O goes to `/mnt/scratch/dasilvaf/iwgc_lncRNA/` (home quota is tight).

---

## 1. Genome + annotation (Phytozome Wm82.a6.v1)

Phytozome requires a free JGI account and accepting the data-usage policy, so this can't be
scripted without your credentials.

**a) Download (web, simplest):** log in at <https://phytozome-next.jgi.doe.gov>, open
*Glycine max Wm82.a6.v1*, and download:
- the **assembly** FASTA — `Gmax_880_v6.0.fa.gz`
- the **annotation** — `Gmax_880_Wm82.a6.v1.gene_exons.gff3.gz`

(exact names are shown in the portal). `scp`/`wget` them onto `amd24`.

**b) Download (scriptable, JGI cookie):**
```bash
curl 'https://signon.jgi.doe.gov/signon/create' \
  --data-urlencode 'login=YOUR_JGI_EMAIL' --data-urlencode 'password=YOUR_PASSWORD' \
  -c /tmp/jgi_cookies > /dev/null
# then download each file with the cookie (copy the file's URL from the portal "Download" link):
curl -b /tmp/jgi_cookies 'https://files.jgi.doe.gov/<path-to>/Gmax_880_v6.0.fa.gz' -o Gmax_genome.fa.gz
curl -b /tmp/jgi_cookies 'https://files.jgi.doe.gov/<path-to>/Gmax_880_Wm82.a6.v1.gene_exons.gff3.gz' -o Gmax_anno.gff3.gz
```

**c) Put them in pipeline format** (`{species}.chromosomes.fa` + `{species}_mRNA.gtf`):
The genome (~1 GB uncompressed) and FASTQ go on **scratch** (home quota). Adjust the source
paths to wherever your Phytozome download landed.
```bash
IN=/mnt/scratch/dasilvaf/iwgc_lncRNA/g.max/01_inputs
SIF=/mnt/scratch/dasilvaf/iwgc_lncRNA/containers/images
export PATH="$PATH:/mnt/ffs24/home/dasilvaf/miniforge3/bin"
mkdir -p "$IN"

gunzip -c .../assembly/Gmax_880_v6.0.fa.gz > "$IN/g.max.chromosomes.fa"

gunzip -c .../annotation/Gmax_880_Wm82.a6.v1.gene_exons.gff3.gz > "$IN/gmax_anno.gff3"
apptainer exec -B /mnt/scratch/dasilvaf/iwgc_lncRNA "$SIF/gffread.sif" \
  gffread "$IN/gmax_anno.gff3" -T -o "$IN/g.max_mRNA.gtf"
rm -f "$IN/gmax_anno.gff3"

# sanity: chromosome names must match between the two files
grep '^>' "$IN/g.max.chromosomes.fa" | head
cut -f1 "$IN/g.max_mRNA.gtf" | grep -v '^#' | sort -u | head
```

> Can't use Phytozome? The closest NCBI substitute is **Glycine_max_v4.0** (`GCF_000004515.6`,
> i.e. Wm82.**a4**, a different/earlier assembly). Prepare it with NCBI `datasets` +
> `scripts/ncbi_datasets_cleanup.py g.max` (as done for Arabidopsis). Document the substitution.

---

## 2. Choose accessions and make the SRA list

The README lists `SRR1174214 SRR1174217 SRR1174218 SRR1174232 …`. Start small:
```bash
printf '%s\n' SRR1174214 SRR1174217 SRR1174218 SRR1174232 > g.max/01_inputs/g.max_sra_IDs.txt
```
(The full 160-accession list is in that same file already if you want everything — but that is a
very large run.)

---

## 3. Download the FASTQ (DEV NODE — compute nodes have no internet)

```bash
bash scripts/download_sra_fastq.sh g.max/01_inputs/g.max_sra_IDs.txt /mnt/scratch/dasilvaf/iwgc_lncRNA/g.max/01_inputs/raw_fastq
```
This pulls `SRRxxxx_1.fastq.gz` + `SRRxxxx_2.fastq.gz` from ENA into `g.max/01_inputs/raw_fastq/`.
~20 GB for the 4 runs — run it in the background if you like:
```bash
nohup bash scripts/download_sra_fastq.sh g.max/01_inputs/g.max_sra_IDs.txt /mnt/scratch/dasilvaf/iwgc_lncRNA/g.max/01_inputs/raw_fastq > dl.log 2>&1 &
tail -f dl.log
```
*(Optional, faithful to the README: download with `prefetch --option-file g.max_sra_IDs.txt` then
`fasterq-dump --split-files` using the `sra_tools.sif` container — also on the dev node.)*

---

## 4. UniProt / DIAMOND DB and containers — already done

- DIAMOND DB: `resources/uniprot/uniprot_out.dmnd` already exists (rebuilt automatically if missing).
- The 15 Singularity images are already built on scratch
  (`/mnt/scratch/dasilvaf/iwgc_lncRNA/containers/images/`). To rebuild:
  `scripts/build_singularity_images.sh`.

---

## 5. Configure

`config/config.gmax.yaml` is ready (outputs/index/containers on scratch). Check two things:
- `strand_specific:` — leave `false`, or set `true` if the libraries are RF strand-specific
  (verify the layout on ENA / SRA; these `SRR1174*` runs are paired-end).
- the accession list / paths point where you put the data.

---

## 6. Dry-run (validates the DAG, no compute)

```bash
conda activate iwgc-lnc-snakemake
snakemake --snakefile workflow/Snakefile --configfile config/config.gmax.yaml \
  --cores 24 --use-singularity --rerun-triggers mtime --dry-run | tail -30
```
You should see jobs for fastp / hisat2_index / hisat2_align / stringtie / … / intersect.

---

## 7. Submit the Slurm job (do NOT run heavy steps on the dev node)

```bash
sbatch scripts/run_pipeline_gmax.sb
```
Defaults: `general-long`, 24 cores, 96 GB, 24 h (soybean genome indexing/alignment is heavy).
Adjust `--time`/`--mem`/`--cpus-per-task` at the top of the script as needed.

---

## 8. Monitor

```bash
squeue -u $USER                              # queue / running state
tail -f iwgc_gmax_<JOBID>.log                # live Snakemake log (in the repo root)
sacct -j <JOBID> --format=JobID,State,Elapsed,MaxRSS,ExitCode   # after it ends
```
Progress lines look like `Finished jobid: N (Rule: ...)` and `X of Y steps (Z%) done`.

---

## 9. Results

```bash
F=/mnt/scratch/dasilvaf/iwgc_lncRNA/g.max/03_outputs/10_final
ls -lah $F
cat $F/Final_lncRNA_results.txt          # high-confidence lncRNA IDs
# also: lncRNA.gtf, lncRNA.fasta, lncRNA_classes.txt, LncRNA_*.txt, Venn_pred_lncRNA.pdf
```

To archive results into the repo (like the Arabidopsis test):
```bash
mkdir -p docs/test_results/gmax_wm82
cp $F/Final_lncRNA_results.txt $F/lncRNA.gtf $F/lncRNA.fasta $F/lncRNA_classes.txt \
   $F/Venn_pred_lncRNA.pdf docs/test_results/gmax_wm82/
```

---

## Troubleshooting

| Symptom | Fix |
|---------|-----|
| `Disk quota exceeded` | you're writing to home; outputs must be on scratch (config already does this). Container images too. VAST quota lags ~30–60 s after `rm`. |
| `prefetch`/download fails in the Slurm job | compute nodes have no internet — download FASTQ on the dev node (step 3), keep `download_sra: false`. |
| `samtools sort: failed to read header from "-"` | HISAT2 must write SAM to stdout (no `-S -`); already handled in the Snakefile. |
| hisat2-build OOM | bump `--mem` in `run_pipeline_gmax.sb` (soybean ~1 Gb). |
| chromosome-name mismatch | genome FASTA headers and GTF column 1 must match (both from Phytozome → they do). |
| resume after a fix | `--rerun-triggers mtime` (already in the sbatch) re-runs only what's missing/changed. |

#!/usr/bin/env Rscript
library(LncFinder)
library(seqinr)   # read input fasta sequences
library(parallel)  # run in parallel
library(e1071)     # for SVM model

# ==========================================
# SET PATHS & CONFIGURATION
# ==========================================
input_fasta   <- "/mnt/scratch/lemasjoh/BASF/IWGC_lncRNA_pipeline/species/03_outputs/04_gffread/species_candidates.fasta"
fold_output   <- "/mnt/scratch/lemasjoh/BASF/IWGC_lncRNA_pipeline/species/03_outputs/08_lncFinder/folded_structures.txt"
output_dir    <- "/mnt/scratch/lemasjoh/BASF/IWGC_lncRNA_pipeline/species/03_outputs/08_lncFinder"
output_csv    <- file.path(output_dir, "species_lncFinder_wheat_model_results.csv")
num_cores     <- 32

cat(">>> [Step 2] Parsing RNAfold structures into memory...\n")

# ==========================================
# PARSE STRUCTURE TEXT FILE
# ==========================================
lines <- readLines(fold_output)

# Identify lines mathematically by position (3 lines per entry in RNAfold output)
header_idx <- seq(1, length(lines), by = 3)
struct_idx <- seq(3, length(lines), by = 3)

# Extract and clean headers
headers <- gsub("^>", "", lines[header_idx])

# Extract dot-bracket notation, dropping the MFE tracking value at the end
raw_structs <- lines[struct_idx]
structures  <- sapply(strsplit(raw_structs, " "), `[`, 1)
names(structures) <- headers

# ==========================================
# LOAD FASTA & ALIGN
# ==========================================
Seqs <- read.fasta(input_fasta, seqtype = "DNA")
fasta_names <- names(Seqs)

if(!all(fasta_names %in% names(structures))) {
  stop("CRITICAL ERROR: Sequence IDs in FASTA and RNAfold output do not match.")
}

# Reorder structure vector to perfectly mirror FASTA arrangement
Structures <- structures[fasta_names]

# ==========================================
# [Step 3] RUN LNCFINDER PREDICTIONS
# ==========================================
cat(">>> [Step 3] Computing LncFinder features and running SVM model...\n")

results <- lnc_finder(
  Seqs,
  Structures       = Structures,
  SS.features      = TRUE,
  format           = "DNA",
  frequencies.file = "wheat",
  svm.model        = "wheat",
  parallel.cores   = num_cores
)

# ==========================================
# SAVE RESULTS
# ==========================================
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

write.csv(results, output_csv, row.names = TRUE)
cat(">>> Success! Results saved to:\n\t", output_csv, "\n")
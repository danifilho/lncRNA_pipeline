#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(LncFinder)
  library(seqinr)
  library(e1071)
})

parse_args <- function(args) {
  opts <- list(threads = 1)
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "-")) {
      stop("Unexpected argument: ", key, call. = FALSE)
    }
    if (i == length(args)) {
      stop("Missing value for argument: ", key, call. = FALSE)
    }
    name <- gsub("-", "_", sub("^--?", "", key))
    opts[[name]] <- args[[i + 1]]
    i <- i + 2
  }
  opts$threads <- as.integer(opts$threads)
  opts
}

opts <- parse_args(commandArgs(trailingOnly = TRUE))

write_empty <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(data.frame(ID = character(), Pred = character()), path, row.names = FALSE)
}

if (is.null(opts$input_fasta) || !file.exists(opts$input_fasta) ||
    file.info(opts$input_fasta)$size == 0) {
  write_empty(opts$output_csv)
  quit(status = 0)
}

lines <- readLines(opts$fold_output, warn = FALSE)
if (!length(lines)) {
  stop("RNAfold output is empty: ", opts$fold_output, call. = FALSE)
}

header_idx <- grep("^>", lines)
if (!length(header_idx)) {
  stop("No FASTA headers were detected in RNAfold output: ", opts$fold_output, call. = FALSE)
}

structure_idx <- header_idx + 2
if (any(structure_idx > length(lines))) {
  stop("RNAfold output does not contain complete three-line records.", call. = FALSE)
}

ids <- sub("^>", "", lines[header_idx])
seq_lines <- lines[header_idx + 1]
struct_lines <- lines[structure_idx]

# Build the LncFinder "SS" data frame: one column per sequence, with three rows
# holding the sequence, the dot-bracket secondary structure, and the minimum
# free energy. This mirrors the layout produced by LncFinder::read_SS() /
# run_RNAfold(), which is what lnc_finder(format = "SS") consumes. The MFE is
# parsed from the trailing "( -12.30)" field exactly as read_SS() does.
build_col <- function(seqstr, structline) {
  n <- nchar(seqstr)
  dotbr <- substr(structline, 1, n)
  mfe <- as.numeric(substr(structline, n + 3, nchar(structline) - 1))
  c(tolower(seqstr), dotbr, mfe)
}

ss_mat <- mapply(build_col, seq_lines, struct_lines)
ss_df <- data.frame(ss_mat, stringsAsFactors = FALSE)
names(ss_df) <- ids

# Some candidate transcripts produce NA features (e.g. no detectable ORF). Left
# in place they silently drop out of predict.svm() and break the internal cbind
# in lnc_finder(). Identify and remove them up front so prediction is stable.
feats <- extract_features(
  ss_df,
  SS.features = TRUE,
  format = "SS",
  frequencies.file = "wheat",
  parallel.cores = opts$threads
)
good_idx <- which(stats::complete.cases(feats))

if (!length(good_idx)) {
  write_empty(opts$output_csv)
  quit(status = 0)
}

ss_df <- ss_df[, good_idx, drop = FALSE]

results <- lnc_finder(
  ss_df,
  SS.features = TRUE,
  format = "SS",
  frequencies.file = "wheat",
  svm.model = "wheat",
  parallel.cores = opts$threads
)

# Emit a tidy ID,Pred table (matches what intersect_ids.R expects: an exact
# "Pred" column with the transcript ID in column 1).
out <- data.frame(
  ID = rownames(results),
  Pred = as.character(results$Pred),
  stringsAsFactors = FALSE
)
dir.create(dirname(opts$output_csv), recursive = TRUE, showWarnings = FALSE)
write.csv(out, opts$output_csv, row.names = FALSE)

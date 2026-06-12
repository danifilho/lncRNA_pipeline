#!/usr/bin/env Rscript

# LncFinder-plant identification, faithful to Plant-LncRNA-pipeline(-v2):
#   frequencies <- make_frequencies(plant training mRNA + lncRNA, SS.features = FALSE)
#   plant       <- readRDS(Plant_model.rda)            # pre-trained e1071 SVM
#   lnc_finder(Seqs, SS.features = FALSE, format = "DNA",
#              frequencies.file = frequencies, svm.model = plant)
# The training FASTAs and the SVM are bundled in the container under
# /opt/Plant-LncRNA-pipline (the upstream repository).

suppressPackageStartupMessages({
  library(LncFinder)
  library(seqinr)
})

args <- commandArgs(trailingOnly = TRUE)
opts <- list(
  threads = 2,
  training_mrna = "/opt/Plant-LncRNA-pipline/example_data/training_mRNA.fasta",
  training_lncrna = "/opt/Plant-LncRNA-pipline/example_data/training_lncRNA.fasta",
  model = "/opt/Plant-LncRNA-pipline/Model/Plant_model.rda"
)
i <- 1
while (i <= length(args)) {
  key <- gsub("-", "_", sub("^--?", "", args[[i]]))
  if (i == length(args)) stop("Missing value for argument: ", args[[i]], call. = FALSE)
  opts[[key]] <- args[[i + 1]]
  i <- i + 2
}
opts$threads <- as.integer(opts$threads)

out <- opts$output_txt
if (is.null(out)) stop("--output-txt is required", call. = FALSE)
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)

# Empty candidate set -> header-only result (keeps the downstream reader happy).
if (is.null(opts$input_fasta) || !file.exists(opts$input_fasta) ||
    file.info(opts$input_fasta)$size == 0) {
  writeLines("Pred\tCoding.Potential", out)
  quit(status = 0)
}

mRNA <- seqinr::read.fasta(file = opts$training_mrna)
lncRNA <- seqinr::read.fasta(file = opts$training_lncrna)

frequencies <- make_frequencies(
  cds.seq = mRNA,
  lncRNA.seq = lncRNA,
  SS.features = FALSE,
  cds.format = "DNA",
  lnc.format = "DNA",
  check.cds = TRUE,
  ignore.illegal = TRUE
)

plant <- readRDS(opts$model)

Seqs <- seqinr::read.fasta(file = opts$input_fasta)
Plant_results <- LncFinder::lnc_finder(
  Seqs,
  SS.features = FALSE,
  format = "DNA",
  frequencies.file = frequencies,
  svm.model = plant,
  parallel.cores = opts$threads
)

write.table(Plant_results, file = out, sep = "\t",
            row.names = TRUE, col.names = TRUE, quote = FALSE)

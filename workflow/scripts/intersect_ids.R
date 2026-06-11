#!/usr/bin/env Rscript

parse_args <- function(args) {
  opts <- list(
    out = "consensus",
    cpat_cutoff = 0.46,
    diamond_evalue = 1e-5
  )
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "-")) {
      stop("Unexpected argument: ", key, call. = FALSE)
    }
    if (grepl("=", key, fixed = TRUE)) {
      parts <- strsplit(sub("^--?", "", key), "=", fixed = TRUE)[[1]]
      name <- gsub("-", "_", parts[[1]])
      value <- parts[[2]]
      i <- i + 1
    } else {
      name <- gsub("-", "_", sub("^--?", "", key))
      if (i == length(args)) {
        stop("Missing value for argument: ", key, call. = FALSE)
      }
      value <- args[[i + 1]]
      i <- i + 2
    }
    opts[[name]] <- value
  }
  opts$cpat_cutoff <- as.numeric(opts$cpat_cutoff)
  opts$diamond_evalue <- as.numeric(opts$diamond_evalue)
  opts
}

read_lines_safe <- function(file) {
  if (is.null(file) || is.na(file) || !nzchar(file) || !file.exists(file) || file.info(file)$size == 0) {
    return(character())
  }
  readLines(file, warn = FALSE)
}

clean_ids <- function(ids) {
  ids <- trimws(as.character(ids))
  ids <- ids[nzchar(ids) & !is.na(ids)]
  unique(ids)
}

read_table_safe <- function(file, sep = "\t", header = TRUE, skip = 0) {
  if (is.null(file) || is.na(file) || !nzchar(file) || !file.exists(file) || file.info(file)$size == 0) {
    return(data.frame())
  }
  tryCatch(
    read.table(
      file,
      sep = sep,
      header = header,
      skip = skip,
      quote = "",
      comment.char = "",
      check.names = FALSE,
      stringsAsFactors = FALSE,
      fill = TRUE
    ),
    error = function(e) data.frame()
  )
}

get_feelnc_ids <- function(file) {
  lines <- read_lines_safe(file)
  if (!length(lines)) return(character())

  if (grepl("\\.gtf$|\\.gff$", file, ignore.case = TRUE)) {
    matches <- regmatches(lines, gregexpr('transcript_id "[^"]+"', lines))
    ids <- unlist(matches, use.names = FALSE)
    ids <- gsub('transcript_id |"', "", ids)
    return(clean_ids(ids))
  }

  clean_ids(lines)
}

get_cpat_ids <- function(file, cutoff) {
  df <- read_table_safe(file, sep = "\t", header = TRUE)
  if (!nrow(df)) return(character())

  names_lower <- tolower(names(df))
  prob_col <- grep("coding.*prob|coding_prob|prob", names_lower)[1]

  # cpat.py writes the sequence ID as an unnamed first column, so read.table
  # (header has one fewer field than the data) puts the IDs into the row names.
  # Use the row names as IDs in that case; otherwise use an explicit seq_ID col.
  if ("seq_id" %in% names_lower) {
    ids <- as.character(df[[which(names_lower == "seq_id")[1]]])
  } else {
    ids <- rownames(df)
  }

  if (is.na(prob_col)) {
    df2 <- read_table_safe(file, sep = "\t", header = FALSE, skip = 1)
    if (!nrow(df2) || ncol(df2) < 6) return(character())
    probs <- suppressWarnings(as.numeric(df2[[6]]))
    return(clean_ids(as.character(df2[[1]])[!is.na(probs) & probs < cutoff]))
  }

  probs <- suppressWarnings(as.numeric(df[[prob_col]]))
  clean_ids(ids[!is.na(probs) & probs < cutoff])
}

get_lncfinder_ids <- function(file) {
  if (is.null(file) || is.na(file) || !nzchar(file) ||
      !file.exists(file) || file.info(file)$size == 0) {
    return(character())
  }
  # LncFinder results are a (quoted) CSV from write.csv(); use read.csv so the
  # surrounding quotes are stripped from both the header and the values. Reading
  # this with read_table_safe(quote = "") would leave literal quotes on the IDs
  # and on "NonCoding", so nothing would ever match below.
  df <- tryCatch(
    read.csv(file, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) data.frame()
  )
  if (!nrow(df)) return(character())

  pred_col <- match("Pred", names(df))
  if (is.na(pred_col)) {
    pred_col <- grep("pred|prediction", names(df), ignore.case = TRUE)[1]
  }
  if (is.na(pred_col)) return(character())

  id_col <- 1
  clean_ids(df[[id_col]][tolower(df[[pred_col]]) %in% c("noncoding", "non-coding", "lncrna", "1")])
}

get_lncboost_ids <- function(file) {
  df <- read_table_safe(file, sep = "\t", header = TRUE)
  if (!nrow(df)) return(character())

  label_col <- match("Predicted_label", names(df))
  if (is.na(label_col)) {
    label_col <- grep("label|pred", names(df), ignore.case = TRUE)[1]
  }
  if (is.na(label_col)) return(character())

  labels <- suppressWarnings(as.numeric(df[[label_col]]))
  clean_ids(df[[1]][!is.na(labels) & labels == 1])
}

get_diamond_hits <- function(file, evalue_cutoff) {
  df <- read_table_safe(file, sep = "\t", header = FALSE)
  if (!nrow(df)) return(character())
  if (ncol(df) >= 11) {
    evalues <- suppressWarnings(as.numeric(df[[11]]))
    return(clean_ids(df[[1]][!is.na(evalues) & evalues <= evalue_cutoff]))
  }
  clean_ids(df[[1]])
}

write_csv_base <- function(df, file) {
  write.csv(df, file = file, row.names = FALSE, quote = TRUE)
}

opts <- parse_args(commandArgs(trailingOnly = TRUE))

feelnc_ids <- get_feelnc_ids(opts$feelnc)
cpat_ids <- get_cpat_ids(opts$cpat, opts$cpat_cutoff)
lncfinder_ids <- get_lncfinder_ids(opts$lncfinder)
boost_ids <- get_lncboost_ids(opts$boost)
protein_hits <- get_diamond_hits(opts$diamond, opts$diamond_evalue)

positive_sets <- list(
  FEELnc = feelnc_ids,
  CPAT = cpat_ids,
  LncFinder = lncfinder_ids,
  Boost = boost_ids
)
provided <- positive_sets[vapply(positive_sets, length, integer(1)) > 0]

if (!length(provided)) {
  final_ids <- character()
  all_ids <- character()
} else {
  final_ids <- Reduce(intersect, provided)
  final_ids <- setdiff(final_ids, protein_hits)
  all_ids <- sort(unique(unlist(provided, use.names = FALSE)))
}

dir.create(dirname(opts$out), recursive = TRUE, showWarnings = FALSE)
writeLines(final_ids, paste0(opts$out, "_final_ids.txt"))

summary_tab <- data.frame(ID = all_ids, stringsAsFactors = FALSE)
summary_tab$Pass_FEELnc <- all_ids %in% feelnc_ids
summary_tab$Pass_CPAT <- all_ids %in% cpat_ids
summary_tab$Pass_LncFinder <- all_ids %in% lncfinder_ids
summary_tab$Pass_Boost <- all_ids %in% boost_ids
summary_tab$Is_Protein <- all_ids %in% protein_hits
summary_tab$Final_Set <- all_ids %in% final_ids
write_csv_base(summary_tab, paste0(opts$out, "_intersection_summary.csv"))

message("Final consensus count: ", length(final_ids))

#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyverse)
  library(optparse)
})

# This custom script accepts outputs from multiple lncRNA prediction tools, extracts candidate IDs, and identifies the consensus set of candidates that pass all filters. It also generates a summary table indicating which candidates passed each tool's criteria.
# This was designed to run with the intersect_ids.sif image.

# 1. DEFINE ARGUMENTS: arguments for each input and output used in the script.
option_list <- list(
  make_option(c("-f", "--feelnc"), type="character", help="FEELnc GTF or TXT list"),
  make_option(c("-c", "--cpat"), type="character", help="CPAT output (TSV)"),
  make_option(c("-b", "--boost"), type="character", help="lncBoost output (TSV)"),
  make_option(c("-l", "--lncfinder"), type="character", help="LncFinder output (CSV)"),
  make_option(c("-d", "--diamond"), type="character", help="DIAMOND output (TSV)"),
  make_option(c("-o", "--out"), type="character", default="consensus", help="Output prefix")
)

opt <- parse_args(OptionParser(option_list=option_list))

# 2. DATA LOADING LOGIC (With Memory Optimization)
get_ids <- function(file, tool_name) {
  if (is.null(file) || !file.exists(file)) {
    cat(tool_name, ": [Skipped]\n")
    return(NULL)
  }
  cat("Loading", tool_name, "... ")
  
  if (tool_name == "FEELnc") {
    if (grepl("\\.gtf$", file, ignore.case = TRUE)) {
      # GTF Parser: Grabs transcript_id "MSTRG.1.1" -> MSTRG.1.1
      lines <- read_lines(file)
      ids <- str_extract(lines, 'transcript_id "[^"]+"') %>%
             str_replace_all('transcript_id |"', '') %>%
             unique()
      ids <- ids[!is.na(ids)]
    } else {
      ids <- read_lines(file)
    }
    
  } else if (tool_name == "CPAT") {
    # Fix header shift: skip header, name columns manually
    df <- read_delim(file, delim = "\t", skip = 1, col_names = FALSE, show_col_types = FALSE)
    ids <- df %>% filter(X6 < 0.46) %>% pull(X1)
    rm(df)
    
  } else if (tool_name == "LncFinder") {
    # Handles "NonCoding" filter and unnamed ID columns
    df <- read_csv(file, show_col_types = FALSE)
    ids <- df %>% filter(Pred == "NonCoding") %>% pull(1)
    rm(df)
    
  } else if (tool_name == "lncBoost") {
    df <- read_delim(file, delim = "\t", show_col_types = FALSE)
    ids <- df %>% filter(Predicted_label == 1) %>% pull(1)
    rm(df)
    
  } else if (tool_name == "DIAMOND") {
    # Negative filter: pulls unique IDs from first column
    df <- read_delim(file, delim = "\t", col_names = FALSE, show_col_types = FALSE)
    ids <- df %>% filter(X11 < 1e-5) %>% pull(1) %>% unique()
    rm(df)
  }
  
  ids <- as.character(trimws(ids))
  cat(length(ids), "candidates identified.\n")
  gc() # Force garbage collection to free RAM
  return(ids)
}

# 3. RUN PIPELINE
feelnc_ids <- get_ids(opt$feelnc, "FEELnc")
cpat_ids   <- get_ids(opt$cpat, "CPAT")
lncf_ids   <- get_ids(opt$lncfinder, "LncFinder")
boost_ids  <- get_ids(opt$boost, "lncBoost")
prot_hits  <- get_ids(opt$diamond, "DIAMOND")

# INTERSECTION
pos_list <- list(FEELnc=feelnc_ids, CPAT=cpat_ids, LncFinder=lncf_ids, Boost=boost_ids)
provided <- pos_list[!sapply(pos_list, is.null)]

if(length(provided) < 1) stop("No input files found!", call.=FALSE)

final_ids <- reduce(provided, intersect)

# Apply DIAMOND filter
if(!is.null(prot_hits)) {
  final_ids <- setdiff(final_ids, prot_hits)
}

# 4. SAVE OUTPUTS
# Save the clean list first (safety in case of OOM on summary)
write_lines(final_ids, paste0(opt$out, "_final_ids.txt"))
cat("\n--- FINAL CONSENSUS COUNT:", length(final_ids), "---\n")

# Build and save Summary Table
cat("Building summary table...\n")
all_ids <- unique(unlist(provided))
summary_tab <- data.frame(ID = all_ids, stringsAsFactors = FALSE)

if(!is.null(feelnc_ids)) summary_tab$Pass_FEELnc  <- all_ids %in% feelnc_ids
if(!is.null(cpat_ids))   summary_tab$Pass_CPAT    <- all_ids %in% cpat_ids
if(!is.null(lncf_ids))   summary_tab$Pass_LncF    <- all_ids %in% lncf_ids
if(!is.null(boost_ids))  summary_tab$Pass_Boost   <- all_ids %in% boost_ids
if(!is.null(prot_hits))  summary_tab$Is_Protein   <- all_ids %in% prot_hits
summary_tab$Final_Set <- all_ids %in% final_ids

write_csv(summary_tab, paste0(opt$out, "_intersection_summary.csv"))
cat("Process complete. Results saved with prefix:", opt$out, "\n")
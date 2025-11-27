#!/usr/bin/env Rscript
library(tidyverse)
source("utils.R")

cosmic_files <- c(
  "COSMIC_v2_SBS_GRCh38.txt",
  "COSMIC_v3.2_SBS_GRCh38.txt"
)

all_results <- list()

# Load input catalogues once
CRC_WES     <- load_crc_data("WES", input_dir, base_dir)
CRC_WGS     <- load_crc_data("WGS", input_dir, base_dir)
CRC_TSO_500 <- load_crc_data("TSO-500", input_dir, base_dir)

for (cf in cosmic_files) {
  message("=== MutationalPatterns with ", cf, " ===")
  
  cosmic_data    <- load_cosmic_signatures(cf, cosmic_dir, base_dir)
  cosmic_version <- gsub("\\.txt$", "", cf)
  
  # Fit
  fit_WES <- fit_crc_with_MutationalPatterns(CRC_WES,     cosmic_data, "WES")
  fit_WGS <- fit_crc_with_MutationalPatterns(CRC_WGS,     cosmic_data, "WGS")
  fit_TSO <- fit_crc_with_MutationalPatterns(CRC_TSO_500, cosmic_data, "TSO-500")
  
  # Cosine
  cos_df <- calculate_cosine_similarity_extended(
    fit_WES, fit_TSO, fit_WGS,
    CRC_WES, CRC_TSO_500, CRC_WGS,
    tool_name      = "MutationalPatterns",
    cosmic_version = cosmic_version
  )
  
  all_results[[cf]] <- cos_df
}

mutp_res <- bind_rows(all_results)

write.table(
  mutp_res,
  file      = "CosineResults_MutationalPatterns_ALL.tsv",
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)

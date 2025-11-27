#!/usr/bin/env Rscript

# Template script:
# - Fits WGS/WES/TSO-500 for each COSMIC file
# - Saves contributions + reconstructed + metadata
# - Outputs a unified CosineResults_*_ALL.tsv

library(tidyverse)
source("utils.R")

## =================== CONFIG ===================

TOOL_NAME <- "MutationalPatterns"  #change name 
fit_fun   <- fit_crc_with_MutationalPatterns  # change per tool

cosmic_files <- c(
  "COSMIC_v2_SBS_GRCh38.txt",
  "COSMIC_v3.2_SBS_GRCh38.txt"
)

# where to save detailed outputs
OUT_BASE <- file.path("signature_results", TOOL_NAME)

# Load input catalogues once
CRC <- list(
  WES      = load_crc_data("WES",     input_dir, base_dir),
  WGS      = load_crc_data("WGS",     input_dir, base_dir),
  `TSO-500` = load_crc_data("TSO-500", input_dir, base_dir)
)

all_cosine <- list()

## =================== MAIN LOOP ===================

for (cf in cosmic_files) {
  message("\n=== ", TOOL_NAME, " with ", cf, " ===")
  
  cosmic_data    <- load_cosmic_signatures(cf, cosmic_dir, base_dir)
  cosmic_version <- gsub("\\.txt$", "", cf)
  
  # output dir for this COSMIC version
  out_dir_version <- file.path(OUT_BASE, cosmic_version)
  dir.create(out_dir_version, recursive = TRUE, showWarnings = FALSE)
  
  # store fits per dataset so we can compute cosine across all three
  fit_list <- list()
  
  for (dt in names(CRC)) {
    message("  -> Fitting ", dt)
    
    fit <- fit_fun(CRC[[dt]], cosmic_data, dt)
    fit_list[[dt]] <- fit
    
    # ---------- Save contributions ----------
    # signatures x samples
    contrib <- as.matrix(fit$contribution)
    write.table(
      contrib,
      file      = file.path(out_dir_version, paste0(dt, "_contribution.tsv")),
      sep       = "\t",
      quote     = FALSE,
      col.names = NA  # keep rownames as first column
    )
    
    # ---------- Save reconstructed catalogues ----------
    # contexts x samples
    recon <- as.matrix(fit$reconstructed)
    write.table(
      recon,
      file      = file.path(out_dir_version, paste0(dt, "_reconstructed.tsv")),
      sep       = "\t",
      quote     = FALSE,
      col.names = NA
    )
    
    # ---------- Save metadata (optional but handy) ----------
    # signatures actually used
    write.table(
      data.frame(Signature = rownames(contrib)),
      file      = file.path(out_dir_version, paste0(dt, "_signatures_used.tsv")),
      sep       = "\t",
      quote     = FALSE,
      row.names = FALSE
    )
    
    # contexts actually used
    write.table(
      data.frame(Context = rownames(recon)),
      file      = file.path(out_dir_version, paste0(dt, "_contexts_used.tsv")),
      sep       = "\t",
      quote     = FALSE,
      row.names = FALSE
    )
  }
  
  # After fitting all three datasets for this COSMIC version:
  # use your extended cosine helper to make tidy rows
  cos_df <- calculate_cosine_similarity_extended(
    fitting_result_WES      = fit_list[["WES"]],
    fitting_result_TSO_500  = fit_list[["TSO-500"]],
    fitting_result_WGS      = fit_list[["WGS"]],
    CRC_WES                 = CRC[["WES"]],
    CRC_TSO_500             = CRC[["TSO-500"]],
    CRC_WGS                 = CRC[["WGS"]],
    tool_name               = TOOL_NAME,
    cosmic_version          = cosmic_version
  )
  
  all_cosine[[cosmic_version]] <- cos_df
}

## =================== SAVE COSINE SUMMARY ===================

cosine_res <- bind_rows(all_cosine)

out_cosine <- paste0("CosineResults_", TOOL_NAME, "_ALL.tsv")
write.table(
  cosine_res,
  file      = out_cosine,
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)

message("\n✅ Saved: ", out_cosine)
message("✅ Detailed contributions + reconstructed profiles saved under: ", OUT_BASE)

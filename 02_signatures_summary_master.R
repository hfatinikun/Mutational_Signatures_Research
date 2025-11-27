#!/usr/bin/env Rscript

# 02_signatures_master.R
#
# Combined script for:
#   - Signature reproducibility across WGS / WES / TSO-500
#   - Signature presence / Jaccard overlap / UpSet plots
#   - COSMIC v2 vs COSMIC v3.2 presence overlap
#
# Replaces:
#   04_signature_reproducibility.R
#   05_signature_overlap_upset.R
#   06_cosmic_version_overlap.R
#   07_plot_cosmic_version_overlap.R


library(tidyverse)
library(pheatmap)


# -------------------- Config --------------------
datasets        <- c("WGS", "WES", "TSO-500")
cosmic_versions <- c("COSMIC_v2_SBS_GRCh38", "COSMIC_v3.2_SBS_GRCh38")
threshold       <- 0.10   # presence threshold (median contribution >= 10%)

base_dirs <- list(
  "MutationalPatterns"    = "signature_results/MutationalPatterns",
  "deconstructSigs"       = "signature_results/deconstructSigs",
  "signature.tools.lib"   = "signature_results/signature.tools.lib",
  "SigLASSO"              = "signature_results/SigLASSO",
  "SigProfilerAssignment" = "signature_results/SPA",
  "MuSiCal"               = "signature_results/Musical"
)

outdir_main       <- "Output/signatures_master"
outdir_repro      <- file.path(outdir_main, "reproducibility")
outdir_overlap    <- file.path(outdir_main, "presence_overlap")
outdir_cosmiccomp <- file.path(outdir_main, "cosmic_version_overlap")

dir.create(outdir_main,       showWarnings = FALSE, recursive = TRUE)
dir.create(outdir_repro,      showWarnings = FALSE, recursive = TRUE)
dir.create(outdir_overlap,    showWarnings = FALSE, recursive = TRUE)
dir.create(outdir_cosmiccomp, showWarnings = FALSE, recursive = TRUE)

# -------------------- Helpers --------------------

# Read normalised contributions if present, else raw and normalise columns
read_contrib_norm_or_raw <- function(dir_version, dataset) {
  p_norm <- file.path(dir_version, paste0(dataset, "_contribution_normalized.tsv"))
  p_raw  <- file.path(dir_version, paste0(dataset, "_contribution.tsv"))
  
  if (file.exists(p_norm)) {
    m <- read.table(p_norm, header = TRUE, sep = "\t",
                    row.names = 1, check.names = FALSE) |> as.matrix()
  } else if (file.exists(p_raw)) {
    m <- read.table(p_raw, header = TRUE, sep = "\t",
                    row.names = 1, check.names = FALSE) |> as.matrix()
    cs <- colSums(m); cs[cs == 0] <- 1
    m <- sweep(m, 2, cs, "/")
  } else {
    return(NULL)
  }
  m
}

# Median exposure per signature (vector)
median_exposure <- function(m) {
  if (is.null(m) || ncol(m) == 0) {
    return(setNames(numeric(nrow(m)), rownames(m)))
  }
  apply(m, 1, median, na.rm = TRUE)
}

# Per-sample correlation across signatures between two datasets
per_sample_corr <- function(matA, matB, pair_name) {
  if (is.null(matA) || is.null(matB)) return(tibble())
  common_sigs <- intersect(rownames(matA), rownames(matB))
  matA <- matA[common_sigs, , drop = FALSE]
  matB <- matB[common_sigs, , drop = FALSE]
  common_samples <- intersect(colnames(matA), colnames(matB))
  if (length(common_samples) == 0) return(tibble())
  
  purrr::map_df(common_samples, function(s) {
    a <- matA[, s, drop = TRUE]
    b <- matB[, s, drop = TRUE]
    tibble(
      Sample        = s,
      corr_pearson  = suppressWarnings(cor(a, b, method = "pearson",  use = "complete.obs")),
      corr_spearman = suppressWarnings(cor(a, b, method = "spearman", use = "complete.obs")),
      pair          = pair_name
    )
  })
}

# Presence from matrix using median >= threshold
presence_from_matrix <- function(m, thr = threshold) {
  if (is.null(m) || ncol(m) == 0) {
    return(tibble(Signature = rownames(m), present = FALSE, value = 0))
  }
  med <- apply(m, 1, median, na.rm = TRUE)
  tibble(
    Signature = names(med),
    present   = as.logical(med >= thr),
    value     = med
  )
}

# Jaccard index
jaccard <- function(a, b) {
  inter <- sum(a & b, na.rm = TRUE)
  union <- sum(a | b, na.rm = TRUE)
  if (union == 0) return(NA_real_)
  inter / union
}

# Proportion of samples with contribution >= thr
prop_samples_ge_thr <- function(m, thr = threshold) {
  if (is.null(m) || ncol(m) == 0) return(rep(0, nrow(m)))
  rowMeans(m >= thr, na.rm = TRUE)
}

# Make UpSet plot via UpSetR if available
make_upset_plot <- function(pres_tbl, out_png) {
  if (!requireNamespace("UpSetR", quietly = TRUE)) {
    message("  (skip UpSet plot: UpSetR not installed)")
    return(invisible(NULL))
  }
  df <- pres_tbl %>%
    select(WGS, WES, `TSO-500`) %>%
    mutate(across(everything(), ~ as.integer(.))) %>%
    as.data.frame()
  rownames(df) <- pres_tbl$Signature
  
  png(out_png, width = 1800, height = 1200, res = 220)
  UpSetR::upset(
    df,
    nsets       = 3,
    sets        = c("WGS", "WES", "TSO-500"),
    nintersects = 10,
    order.by    = "freq"
  )
  dev.off()
}

# -------------------- Global accumulators --------------------

# For dataset-level median-exposure corr heatmap (averaged across tools & COSMIC)
dataset_corr_sum <- matrix(
  0,
  nrow = length(datasets),
  ncol = length(datasets),
  dimnames = list(datasets, datasets)
)

dataset_corr_count <- matrix(
  0,
  nrow = length(datasets),
  ncol = length(datasets),
  dimnames = list(datasets, datasets)
)

# For per-sample corr summary
all_samplecorr_rows   <- list()

# For presence-based COSMIC v2 vs v3.2 overlap (will reuse presence tables)
presence_files <- list()  # store presence file paths keyed by Tool & COSMIC

# -------------------- Main loop: Tool × COSMIC --------------------

for (tool in names(base_dirs)) {
  for (cosmic in cosmic_versions) {
    message("\n=== ", tool, " | ", cosmic, " ===")
    
    dir_version <- file.path(base_dirs[[tool]], cosmic)
    if (!dir.exists(dir_version)) {
      message("  (skip: not found) ", dir_version)
      next
    }
    
    # ---- Load matrices for the three datasets ----
    mats <- setNames(vector("list", length(datasets)), datasets)
    for (dt in datasets) {
      mats[[dt]] <- read_contrib_norm_or_raw(dir_version, dt)
    }
    if (any(sapply(mats, is.null))) {
      message("  (skip: missing at least one dataset)")
      next
    }
    
    # Align signatures across datasets
    common_sigs <- Reduce(intersect, lapply(mats, rownames))
    mats <- lapply(mats, function(m) m[common_sigs, , drop = FALSE])
    
    # ============================================================
    # A) REPRODUCIBILITY: dataset-level + per-sample
    # ============================================================
    
    # Median exposure per dataset
    med_list <- lapply(mats, median_exposure)
    med_df <- tibble(Signature = common_sigs)
    for (dt in datasets) {
      med_df[[dt]] <- as.numeric(med_list[[dt]])
    }
    
    # 3×3 dataset correlation matrix (median exposures)
    med_mat <- med_df %>% select(all_of(datasets))
    ds_corr <- suppressWarnings(
      cor(med_mat, method = "spearman", use = "pairwise.complete.obs")
    )
    
    # Accumulate for global average heatmap
    valid <- !is.na(ds_corr)
    dataset_corr_sum[valid]   <- dataset_corr_sum[valid]   + ds_corr[valid]
    dataset_corr_count[valid] <- dataset_corr_count[valid] + 1
    
    # Save per-tool/cosmic dataset corr as a small table
    ds_corr_long <- as.data.frame(as.table(ds_corr))
    colnames(ds_corr_long) <- c("Dataset1", "Dataset2", "Spearman")
    
    write.table(
      ds_corr_long,
      file.path(
        outdir_repro,
        paste0("DatasetCorr_", tool, "_", cosmic, ".tsv")
      ),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
    
    # Per-sample correlations across signatures
    ps_WGS_WES <- per_sample_corr(mats$WGS, mats$WES,       "WGS_vs_WES")
    ps_WES_TSO <- per_sample_corr(mats$WES, mats$`TSO-500`, "WES_vs_TSO-500")
    ps_WGS_TSO <- per_sample_corr(mats$WGS, mats$`TSO-500`, "WGS_vs_TSO-500")
    ps_all <- bind_rows(ps_WGS_WES, ps_WES_TSO, ps_WGS_TSO) %>%
      mutate(Tool = tool, COSMIC_Version = cosmic)
    
    if (nrow(ps_all) > 0) {
      write.table(
        ps_all,
        file.path(
          outdir_repro,
          paste0("PerSampleCorr_", tool, "_", cosmic, ".tsv")
        ),
        sep = "\t", quote = FALSE, row.names = FALSE
      )
      all_samplecorr_rows[[paste(tool, cosmic)]] <- ps_all
    }
    
    # ============================================================
    # B) PRESENCE / JACCARD / UPSET / HEATMAP
    # ============================================================
    
    # Presence table (median >= threshold)
    pres_list <- lapply(datasets, function(dt) {
      tmp <- presence_from_matrix(mats[[dt]], thr = threshold) %>%
        select(Signature, present)
      colnames(tmp)[2] <- dt
      tmp
    })
    
    pres_tbl <- purrr::reduce(pres_list, full_join, by = "Signature") %>%
      replace_na(list(WGS = FALSE, WES = FALSE, `TSO-500` = FALSE))
    
    pres_out <- file.path(
      outdir_overlap,
      paste0(
        "Presence_", tool, "_", cosmic,
        "_median_thr", threshold * 100, ".tsv"
      )
    )
    write.table(
      pres_tbl,
      pres_out,
      sep = "\t", quote = FALSE, row.names = FALSE
    )
    
    # Remember presence files for COSMIC v2 vs v3.2 comparison later
    presence_files[[paste(tool, cosmic)]] <- pres_out
    
    # Jaccard indices between datasets
    J <- tibble(
      pair    = c("WGS_vs_WES", "WES_vs_TSO-500", "WGS_vs_TSO-500"),
      jaccard = c(
        jaccard(pres_tbl$WGS, pres_tbl$WES),
        jaccard(pres_tbl$WES, pres_tbl$`TSO-500`),
        jaccard(pres_tbl$WGS, pres_tbl$`TSO-500`)
      )
    )
    
    jac_out <- file.path(
      outdir_overlap,
      paste0(
        "Jaccard_", tool, "_", cosmic,
        "_median_thr", threshold * 100, ".tsv"
      )
    )
    write.table(
      J,
      jac_out,
      sep = "\t", quote = FALSE, row.names = FALSE
    )
    
    # UpSet plot
    make_upset_plot(
      pres_tbl,
      file.path(
        outdir_overlap,
        paste0(
          "UpSet_", tool, "_", cosmic,
          "_median_thr", threshold * 100, ".png"
        )
      )
    )
    
    # Proportion-of-samples heatmap
    prop_df <- tibble(
      Signature = common_sigs,
      WGS       = prop_samples_ge_thr(mats$WGS,      thr = threshold),
      WES       = prop_samples_ge_thr(mats$WES,      thr = threshold),
      `TSO-500` = prop_samples_ge_thr(mats$`TSO-500`, thr = threshold)
    ) %>%
      column_to_rownames("Signature")
    
    pheatmap(
      as.matrix(prop_df),
      color         = colorRampPalette(c("white", "grey80", "black"))(100),
      cluster_rows  = TRUE,
      cluster_cols  = FALSE,
      breaks        = seq(0, 1, length.out = 101),
      main          = paste0(
        tool, " | ", cosmic,
        " | proportion of samples ≥ ", threshold * 100, "%"
      ),
      filename      = file.path(
        outdir_overlap,
        paste0(
          "Heatmap_PropSamples_", tool, "_", cosmic,
          "_thr", threshold * 100, ".png"
        )
      ),
      width         = 6,
      height        = 10
    )
  }
}

# ============================================================
# C) GLOBAL SUMMARY: DATASET-LEVEL HEATMAP & PER-SAMPLE SUMMARY
# ============================================================

# Global dataset correlation heatmap (averaged across tools & COSMIC)
if (sum(dataset_corr_count) > 0) {
  avg_corr <- dataset_corr_sum / dataset_corr_count
  
  pheatmap(
    avg_corr,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    main         = "Median exposure corr (datasets)\nAveraged across tools & COSMIC versions",
    filename     = file.path(outdir_repro, "Heatmap_DatasetMedianCorr_AllTools_AllCOSMIC.png"),
    width        = 5,
    height       = 4
  )
}

# Per-sample correlation summary across tools & COSMIC
if (length(all_samplecorr_rows) > 0) {
  per_sample_summary <- bind_rows(all_samplecorr_rows) %>%
    group_by(Tool, COSMIC_Version, pair) %>%
    summarise(
      n              = n(),
      median_spearman = median(corr_spearman, na.rm = TRUE),
      iqr_spearman    = IQR(corr_spearman, na.rm = TRUE),
      .groups         = "drop"
    )
  
  write.table(
    per_sample_summary,
    file.path(outdir_repro, "PerSampleCorr_Summary_byToolCOSMIC.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}

# ============================================================
# D) COSMIC v2 vs v3.2 OVERLAP (PRESENCE-BASED) + PLOT
# ============================================================

tools <- names(base_dirs)
cosmic_v2  <- "COSMIC_v2_SBS_GRCh38"
cosmic_v32 <- "COSMIC_v3.2_SBS_GRCh38"

results <- list()

for (tool in tools) {
  message("\n[ COSMIC overlap ] ", tool)
  
  key_v2  <- paste(tool, cosmic_v2)
  key_v32 <- paste(tool, cosmic_v32)
  
  if (! (key_v2 %in% names(presence_files) && key_v32 %in% names(presence_files)) ) {
    message("  (skip: missing presence files for both COSMIC versions)")
    next
  }
  
  pres_v2  <- read.table(
    presence_files[[key_v2]],
    header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE
  )
  pres_v32 <- read.table(
    presence_files[[key_v32]],
    header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE
  )
  
  all_sigs <- union(pres_v2$Signature, pres_v32$Signature)
  pres_v2_full  <- pres_v2  %>% right_join(tibble(Signature = all_sigs), by = "Signature")
  pres_v32_full <- pres_v32 %>% right_join(tibble(Signature = all_sigs), by = "Signature")
  
  pres_v2_full[datasets]  <- lapply(pres_v2_full[datasets],  function(x) ifelse(is.na(x), FALSE, x))
  pres_v32_full[datasets] <- lapply(pres_v32_full[datasets], function(x) ifelse(is.na(x), FALSE, x))
  
  for (dt in datasets) {
    v2_vec  <- as.logical(pres_v2_full[[dt]])
    v32_vec <- as.logical(pres_v32_full[[dt]])
    
    n_v2       <- sum(v2_vec,  na.rm = TRUE)
    n_v32      <- sum(v32_vec, na.rm = TRUE)
    both_vec   <- v2_vec & v32_vec
    n_both     <- sum(both_vec, na.rm = TRUE)
    n_v2_only  <- sum(v2_vec  & !v32_vec, na.rm = TRUE)
    n_v32_only <- sum(!v2_vec &  v32_vec, na.rm = TRUE)
    jac        <- jaccard(v2_vec, v32_vec)
    
    results[[length(results) + 1]] <- tibble(
      Tool          = tool,
      Data_Type     = dt,
      COSMIC_v2     = cosmic_v2,
      COSMIC_v3.2   = cosmic_v32,
      threshold     = threshold,
      n_v2_present  = n_v2,
      n_v32_present = n_v32,
      n_both        = n_both,
      n_v2_only     = n_v2_only,
      n_v32_only    = n_v32_only,
      jaccard       = jac
    )
  }
}

if (length(results) > 0) {
  summary_tbl <- bind_rows(results)
  
  out_file <- file.path(outdir_cosmiccomp, "COSMIC_v2_vs_v3.2_signature_overlap_summary.tsv")
  write.table(
    summary_tbl,
    out_file,
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  # Plot Jaccard overlap
  summary_tbl$Tool      <- factor(summary_tbl$Tool)
  summary_tbl$Data_Type <- factor(summary_tbl$Data_Type,
                                  levels = c("WGS", "WES", "TSO-500"))
  
  p_cos <- ggplot(summary_tbl,
                  aes(x = jaccard,
                      y = Tool,
                      colour = Data_Type)) +
    geom_point(size = 3) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    theme_bw() +
    labs(
      title = "Overlap of signatures: COSMIC v2 vs v3.2",
      x     = "Jaccard index (presence v2 vs v3.2)",
      y     = "Tool"
    )
  
  ggsave(
    file.path(outdir_cosmiccomp, "Figure_COSMIC_v2_vs_v3.2_Jaccard.png"),
    p_cos,
    width = 7, height = 4, dpi = 300
  )
}

message("\n✅ 02_signatures_master.R completed.")
message("  Reproducibility outputs:      ", outdir_repro)
message("  Presence/overlap outputs:    ", outdir_overlap)
message("  COSMIC v2 vs v3.2 outputs:   ", outdir_cosmiccomp)

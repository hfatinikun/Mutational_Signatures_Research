#!/usr/bin/env Rscript

# 03_extended_analyses_master.R
#
# Extended analyses on top of fitted signature contributions:
#   1) Tool-to-tool reproducibility (exposure correlations)
#   2) Fragile signatures (unstable across tools)
#   3) Signature co-occurrence (signature–signature correlations)
#   4) Sample clustering (PCA of signature exposures)
#
# Uses: signature_results/<Tool>/<COSMIC>/<DataType>_contribution[_normalized].tsv

library(tidyverse)
library(pheatmap)


# -------------------- Config --------------------

datasets        <- c("WGS", "WES", "TSO-500")
cosmic_versions <- c("COSMIC_v2_SBS_GRCh38", "COSMIC_v3.2_SBS_GRCh38")

tools <- c(
  "MutationalPatterns",
  "deconstructSigs",
  "signature.tools.lib",
  "SigLASSO",
  "SPA",       # adjust to match your actual folder name if needed
  "Musical"    # adjust if different
)

# Map tool name -> directory (adapt if your folder names differ)
base_dirs <- list(
  "MutationalPatterns"  = "signature_results/MutationalPatterns",
  "deconstructSigs"     = "signature_results/deconstructSigs",
  "signature.tools.lib" = "signature_results/signature.tools.lib",
  "SigLASSO"            = "signature_results/SigLASSO",
  "SPA"                 = "signature_results/SPA",
  "Musical"             = "signature_results/Musical"
)

outdir_main          <- "Output/extended_analyses"
outdir_tool2tool     <- file.path(outdir_main, "tool_to_tool")
outdir_fragile       <- file.path(outdir_main, "fragile_signatures")
outdir_cooccur       <- file.path(outdir_main, "cooccurrence")
outdir_clustering    <- file.path(outdir_main, "clustering")

dir.create(outdir_main,      showWarnings = FALSE, recursive = TRUE)
dir.create(outdir_tool2tool, showWarnings = FALSE, recursive = TRUE)
dir.create(outdir_fragile,   showWarnings = FALSE, recursive = TRUE)
dir.create(outdir_cooccur,   showWarnings = FALSE, recursive = TRUE)
dir.create(outdir_clustering,showWarnings = FALSE, recursive = TRUE)

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

# Compute Spearman corr between two numeric vectors, safe
safe_spearman <- function(x, y) {
  suppressWarnings(cor(x, y, method = "spearman", use = "complete.obs"))
}

# -------------------- 1) TOOL-TO-TOOL REPRODUCIBILITY --------------------
# For each dataset × COSMIC × signature, compute Spearman correlation of exposures
# between every pair of tools.

message("\n=== 1) Tool-to-tool reproducibility ===")

toolpair_results <- list()

for (dt in datasets) {
  for (cosmic in cosmic_versions) {
    message("  Dataset = ", dt, " | COSMIC = ", cosmic)
    
    # Read contributions for all tools that have this combo
    mats <- list()
    for (tool in tools) {
      dir_version <- file.path(base_dirs[[tool]], cosmic)
      if (!dir.exists(dir_version)) next
      m <- read_contrib_norm_or_raw(dir_version, dt)
      if (!is.null(m) && ncol(m) > 0) {
        mats[[tool]] <- m
      }
    }
    
    if (length(mats) < 2) {
      message("    (skip: <2 tools available)")
      next
    }
    
    # Align signatures across all tools
    common_sigs <- Reduce(intersect, lapply(mats, rownames))
    if (length(common_sigs) == 0) {
      message("    (skip: no common signatures across tools)")
      next
    }
    mats <- lapply(mats, function(m) m[common_sigs, , drop = FALSE])
    
    # Align samples per tool pair at correlation step
    tool_names <- names(mats)
    if (length(tool_names) < 2) next
    
    # All unordered tool pairs
    tool_pairs <- t(combn(tool_names, 2))
    
    for (k in seq_len(nrow(tool_pairs))) {
      t1 <- tool_pairs[k, 1]
      t2 <- tool_pairs[k, 2]
      m1 <- mats[[t1]]
      m2 <- mats[[t2]]
      
      # Align samples for this pair
      common_samples <- intersect(colnames(m1), colnames(m2))
      if (length(common_samples) < 2) next   # need at least 2 samples
      
      m1_sub <- m1[, common_samples, drop = FALSE]
      m2_sub <- m2[, common_samples, drop = FALSE]
      
      # For each signature, correlation across samples
      sig_corr <- purrr::map_dbl(common_sigs, function(sig) {
        safe_spearman(m1_sub[sig, ], m2_sub[sig, ])
      })
      
      df_pair <- tibble(
        Data_Type      = dt,
        COSMIC_Version = cosmic,
        Tool1          = t1,
        Tool2          = t2,
        Signature      = common_sigs,
        Spearman       = sig_corr
      )
      
      toolpair_results[[length(toolpair_results) + 1]] <- df_pair
    }
  }
}

if (length(toolpair_results) > 0) {
  toolpair_all <- bind_rows(toolpair_results)
  
  # Save full per-signature table
  write.table(
    toolpair_all,
    file.path(outdir_tool2tool, "ToolPair_SignatureCorr_All.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  # Summary: median Spearman per tool pair / dataset / COSMIC
  toolpair_summary <- toolpair_all %>%
    group_by(Data_Type, COSMIC_Version, Tool1, Tool2) %>%
    summarise(
      n_signatures   = sum(!is.na(Spearman)),
      median_spearman = median(Spearman, na.rm = TRUE),
      iqr_spearman    = IQR(Spearman, na.rm = TRUE),
      .groups         = "drop"
    )
  
  write.table(
    toolpair_summary,
    file.path(outdir_tool2tool, "ToolPair_SignatureCorr_Summary.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  # Quick global plot: median Spearman per tool pair (faceted by dataset)
  p_tp <- ggplot(toolpair_summary,
                 aes(x = median_spearman,
                     y = interaction(Tool1, Tool2),
                     colour = Data_Type)) +
    geom_point(size = 2) +
    facet_wrap(~ COSMIC_Version) +
    theme_bw() +
    labs(
      title = "Tool-to-tool reproducibility (signature exposures)",
      x     = "Median Spearman (per signature, across samples)",
      y     = "Tool pair"
    )
  
  ggsave(
    file.path(outdir_tool2tool, "Figure_ToolPair_SignatureCorr_Summary.png"),
    p_tp, width = 9, height = 5, dpi = 300
  )
  
} else {
  message("⚠ No tool-to-tool correlations computed.")
}

# -------------------- 2) FRAGILE SIGNATURES --------------------
# Signatures with low median correlation across tool pairs (unstable across tools)

message("\n=== 2) Fragile signatures ===")

if (exists("toolpair_all")) {
  fragile_tbl <- toolpair_all %>%
    group_by(Signature) %>%
    summarise(
      n_values        = sum(!is.na(Spearman)),
      global_median   = median(Spearman, na.rm = TRUE),
      global_IQR      = IQR(Spearman, na.rm = TRUE),
      .groups         = "drop"
    ) %>%
    arrange(global_median)
  
  write.table(
    fragile_tbl,
    file.path(outdir_fragile, "Signature_ToolStability_Global.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  # Define "fragile" as global_median < 0.7 (you can change this)
  fragile_cutoff <- 0.7
  fragile_only <- fragile_tbl %>%
    filter(!is.na(global_median), global_median < fragile_cutoff)
  
  write.table(
    fragile_only,
    file.path(outdir_fragile, paste0("FragileSignatures_median_lt_", fragile_cutoff, ".tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
} else {
  message("⚠ No toolpair_all found; skip fragile signatures section.")
}

# -------------------- 3) SIGNATURE CO-OCCURRENCE --------------------
# Correlations between signatures across samples (same tool & COSMIC).
# We pick a "reference tool" and COSMIC v3.2 for cleaner story.

message("\n=== 3) Signature co-occurrence ===")

ref_tool   <- "MutationalPatterns"       # adjust if you want a different base tool
ref_cosmic <- "COSMIC_v3.2_SBS_GRCh38"   # focus on v3.2

cooccur_results <- list()

for (dt in datasets) {
  dir_version <- file.path(base_dirs[[ref_tool]], ref_cosmic)
  if (!dir.exists(dir_version)) {
    message("  (skip co-occurrence: dir not found) ", dir_version)
    next
  }
  m <- read_contrib_norm_or_raw(dir_version, dt)
  if (is.null(m) || ncol(m) < 3) {
    message("  (skip co-occurrence: not enough samples) for ", dt)
    next
  }
  
  # Filter out signatures that are almost never used
  row_means <- rowMeans(m, na.rm = TRUE)
  keep_sigs <- names(row_means)[row_means > 0.01]  # threshold can be tuned
  if (length(keep_sigs) < 2) {
    message("  (skip co-occurrence: <2 active signatures) for ", dt)
    next
  }
  m_sub <- m[keep_sigs, , drop = FALSE]
  
  # Corr across samples: signatures x signatures (Spearman)
  sig_cor <- suppressWarnings(
    cor(t(m_sub), method = "spearman", use = "pairwise.complete.obs")
  )
  
  # Save matrix and heatmap
  out_mat <- file.path(outdir_cooccur,
                       paste0("Cooccurrence_SignatureCorr_", ref_tool, "_", ref_cosmic, "_", dt, ".tsv"))
  write.table(
    sig_cor,
    out_mat,
    sep = "\t", quote = FALSE, row.names = TRUE
  )
  
  pheatmap(
    sig_cor,
    main     = paste0("Signature co-occurrence: ", ref_tool, " | ", ref_cosmic, " | ", dt),
    filename = file.path(outdir_cooccur,
                         paste0("Heatmap_Cooccurrence_", ref_tool, "_", ref_cosmic, "_", dt, ".png")),
    width    = 8,
    height   = 8
  )
  
  cooccur_results[[dt]] <- sig_cor
}

# -------------------- 4) SAMPLE CLUSTERING (PCA) --------------------
# Combine samples across datasets for a reference tool/COSMIC and run PCA.

message("\n=== 4) Sample clustering (PCA on exposures) ===")

pca_tool   <- "MutationalPatterns"
pca_cosmic <- "COSMIC_v3.2_SBS_GRCh38"

# Collect contributions for each dataset & tag samples with Data_Type
expo_list <- list()

for (dt in datasets) {
  dir_version <- file.path(base_dirs[[pca_tool]], pca_cosmic)
  if (!dir.exists(dir_version)) next
  m <- read_contrib_norm_or_raw(dir_version, dt)
  if (is.null(m) || ncol(m) == 0) next
  
  # Transpose: samples × signatures
  df <- as.data.frame(t(m))
  df$Sample    <- rownames(df)
  df$Data_Type <- dt
  expo_list[[dt]] <- df
}

if (length(expo_list) > 0) {
  expo_all <- bind_rows(expo_list)
  
  # Keep numeric signature columns only
  sig_cols <- setdiff(colnames(expo_all), c("Sample", "Data_Type"))
  X <- expo_all[, sig_cols, drop = FALSE]
  
  # Remove columns with zero variance
  keep_cols <- apply(X, 2, function(v) sd(v, na.rm = TRUE) > 0)
  X <- X[, keep_cols, drop = FALSE]
  
  if (ncol(X) >= 2) {
    # PCA on scaled exposures
    pca_res <- prcomp(X, center = TRUE, scale. = TRUE)
    
    scores <- as.data.frame(pca_res$x[, 1:2, drop = FALSE])
    scores$Sample    <- expo_all$Sample
    scores$Data_Type <- expo_all$Data_Type
    
    write.table(
      scores,
      file.path(outdir_clustering,
                paste0("PCA_scores_", pca_tool, "_", pca_cosmic, ".tsv")),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
    
    var_expl <- (pca_res$sdev^2) / sum(pca_res$sdev^2)
    pc1_lab  <- paste0("PC1 (", round(var_expl[1] * 100, 1), "%)")
    pc2_lab  <- paste0("PC2 (", round(var_expl[2] * 100, 1), "%)")
    
    p_pca <- ggplot(scores, aes(x = PC1, y = PC2, colour = Data_Type)) +
      geom_point(alpha = 0.7, size = 2) +
      theme_bw() +
      labs(
        title = paste0("PCA of signature exposures: ", pca_tool, " | ", pca_cosmic),
        x     = pc1_lab,
        y     = pc2_lab
      )
    
    ggsave(
      file.path(outdir_clustering,
                paste0("PCA_", pca_tool, "_", pca_cosmic, ".png")),
      p_pca, width = 7, height = 5, dpi = 300
    )
  } else {
    message("  (skip PCA: not enough varying signatures)")
  }
} else {
  message("  (skip PCA: no exposures available)")
}

message("\n✅ 03_extended_analyses_master.R completed.")
message("  Tool-to-tool outputs:        ", outdir_tool2tool)
message("  Fragile signature outputs:   ", outdir_fragile)
message("  Co-occurrence outputs:       ", outdir_cooccur)
message("  Clustering outputs:          ", outdir_clustering)

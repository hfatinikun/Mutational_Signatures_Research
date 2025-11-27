#!/usr/bin/env Rscript

## ============================================================
## 04_master_additional_analyses.R
##
## Performs:
##  1) Minimum mutation threshold analysis (subsampling)
##  2) ΔMMR / ΔPOLE (reference-only) from signature contributions
##  3) WGS vs WES vs TSO coding vs whole-genome comparison
##  4) Flat signature inflation test (via WGS vs WES)
##
## Uses tools: MutationalPatterns, SigLASSO
## Uses COSMIC: v2 and v3.2 (GRCh38)
##
## Output root: Output/master_analysis/
## ============================================================


library(tidyverse)
library(readr)
library(tibble)


## ------------------------------------------------------------
## Setup
## ------------------------------------------------------------

base_dir   <- getwd()
input_dir  <- file.path(base_dir, "Input files")
sig_dir    <- file.path(base_dir, "signature_results")
out_root   <- file.path(base_dir, "Output", "master_analysis")

if (!dir.exists(out_root)) dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

cosmic_versions <- c("COSMIC_v2_SBS_GRCh38", "COSMIC_v3.2_SBS_GRCh38")
tools_of_interest <- c("MutationalPatterns", "SigLASSO")
datasets_preclinical <- c("WES", "WGS", "TSO-500")  # as in your signature_results structure

## Helper for safe message
msg <- function(...) cat("[MASTER]", ..., "\n")

## ------------------------------------------------------------
## Utility functions
## ------------------------------------------------------------

cosine_similarity <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  num <- sum(x * y)
  den <- sqrt(sum(x^2)) * sqrt(sum(y^2))
  if (den == 0) return(NA_real_)
  num / den
}

## Read catalogue matrix (contexts x samples)
read_catalogue <- function(path) {
  read.table(path, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
}

## Subsample a profile (vector of counts) at given fraction
subsample_profile <- function(counts, frac) {
  counts <- as.numeric(counts)
  total <- sum(counts)
  if (total <= 0) return(rep(0, length(counts)))
  
  n_draw <- max(1, floor(total * frac))
  # Expand into vector of mutation indices
  idx_vec <- rep(seq_along(counts), counts)
  drawn <- sample(idx_vec, size = n_draw, replace = FALSE)
  tabulate(drawn, nbins = length(counts))
}

## Load contribution matrix for a given tool / cosmic / dataset
## Assumes signatures as rows, samples as columns
load_contributions_long <- function(tool, cosmic_version, dataset) {
  # Prefer normalized file if present, else raw and normalize
  contrib_base <- file.path(sig_dir, tool, cosmic_version)
  norm_path <- file.path(contrib_base, paste0(dataset, "_contribution_normalized.tsv"))
  raw_path  <- file.path(contrib_base, paste0(dataset, "_contribution.tsv"))
  
  if (file.exists(norm_path)) {
    df <- read.table(norm_path, header = TRUE, sep = "\t",
                     check.names = FALSE, row.names = 1)
    long <- df %>%
      rownames_to_column("signature") %>%
      pivot_longer(-signature, names_to = "sample", values_to = "fraction")
  } else if (file.exists(raw_path)) {
    df <- read.table(raw_path, header = TRUE, sep = "\t",
                     check.names = FALSE, row.names = 1)
    long <- df %>%
      rownames_to_column("signature") %>%
      pivot_longer(-signature, names_to = "sample", values_to = "raw_contribution") %>%
      group_by(sample) %>%
      mutate(fraction = ifelse(sum(raw_contribution) > 0,
                               raw_contribution / sum(raw_contribution),
                               0)) %>%
      ungroup()
  } else {
    msg("WARNING: No contribution file found for",
        tool, cosmic_version, dataset, " — skipping.")
    return(NULL)
  }
  
  long %>%
    mutate(tool = tool,
           cosmic_version = cosmic_version,
           dataset = dataset) %>%
    select(tool, cosmic_version, dataset, sample, signature, fraction)
}

## ------------------------------------------------------------
## 1. Minimum mutation threshold analysis
## ------------------------------------------------------------

msg("Running minimum mutation threshold analysis...")

min_mut_outfile <- file.path(out_root, "min_mutation_threshold.tsv")

# Input catalogue files
catalogue_files <- tibble(
  file = c("Preclinical_Dataset_WES.txt",
           "Preclinical_Dataset_WGS.txt",
           "Preclinical_Dataset_TSO-500.txt",
           "Clinical_Dataset_WES.txt"),
  dataset_label = c("Preclinical_WES",
                    "Preclinical_WGS",
                    "Preclinical_TSO-500",
                    "Clinical_WES")
)

subsample_fracs <- seq(0.05, 0.95, by = 0.05)
n_rep <- 5  # number of replicates per fraction

threshold_results <- list()

for (i in seq_len(nrow(catalogue_files))) {
  file_name <- catalogue_files$file[i]
  dataset   <- catalogue_files$dataset_label[i]
  path <- file.path(input_dir, file_name)
  
  if (!file.exists(path)) {
    msg("WARNING: Catalogue file not found:", path, " — skipping.")
    next
  }
  
  msg("  Processing catalogue:", file_name, "(", dataset, ")")
  
  mat <- read_catalogue(path)  # contexts x samples
  
  for (sample_name in colnames(mat)) {
    full_profile <- mat[, sample_name]
    full_total   <- sum(full_profile)
    
    if (full_total <= 0) {
      msg("   - Sample", sample_name, "has zero mutations, skipping.")
      next
    }
    
    for (f in subsample_fracs) {
      for (r in seq_len(n_rep)) {
        subs <- subsample_profile(full_profile, f)
        cs   <- cosine_similarity(subs, full_profile)
        threshold_results[[length(threshold_results) + 1]] <- tibble(
          dataset        = dataset,
          sample         = sample_name,
          fraction       = f,
          replicate      = r,
          n_mut_full     = full_total,
          n_mut_subsample = sum(subs),
          cosine         = cs
        )
      }
    }
  }
}

if (length(threshold_results) > 0) {
  threshold_df <- bind_rows(threshold_results)
  write_tsv(threshold_df, min_mut_outfile)
  msg("  -> Written subsampling results to:", min_mut_outfile)
  
  # Summary: median cosine per dataset x mutation bin
  threshold_summary <- threshold_df %>%
    group_by(dataset, fraction) %>%
    summarise(
      median_cosine = median(cosine, na.rm = TRUE),
      median_n_mut  = median(n_mut_subsample, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Plot cosine vs number of mutations (log10 x-axis optional)
  p_thr <- ggplot(threshold_summary,
                  aes(x = median_n_mut, y = median_cosine,
                      colour = dataset)) +
    geom_line() +
    geom_point() +
    theme_bw() +
    labs(x = "Median number of mutations (subsampled)",
         y = "Median cosine similarity\n(subsample vs full profile)",
         title = "Minimum mutation threshold analysis")
  
  ggsave(file.path(out_root, "Figure_MinMut_vs_Cosine.png"),
         p_thr, width = 7, height = 5, dpi = 300)
  
  msg("  -> Plot saved: Figure_MinMut_vs_Cosine.png")
} else {
  msg("  No subsampling results generated.")
}

## ------------------------------------------------------------
## 2. Load contributions for MP & SigLASSO (all datasets)
## ------------------------------------------------------------

msg("Loading contribution matrices (MP & SigLASSO)...")

all_contrib <- list()

for (tool in tools_of_interest) {
  for (cv in cosmic_versions) {
    for (ds in datasets_preclinical) {
      df_long <- load_contributions_long(tool, cv, ds)
      if (!is.null(df_long)) {
        all_contrib[[length(all_contrib) + 1]] <- df_long
      }
    }
  }
}

if (length(all_contrib) == 0) {
  stop("No contribution files found for MutationalPatterns/SigLASSO — aborting.")
}

contrib_df <- bind_rows(all_contrib)

# Ensure factors
contrib_df <- contrib_df %>%
  mutate(
    tool = factor(tool, levels = tools_of_interest),
    dataset = factor(dataset, levels = datasets_preclinical)
  )

## ------------------------------------------------------------
## 3. ΔMMR / ΔPOLE (reference-only version)
## ------------------------------------------------------------

msg("Computing ΔMMR / ΔPOLE (reference-only)...")

# Signature definitions (you can tweak these if needed)
mmr_sigs_v2   <- c("SBS6", "SBS15", "SBS20", "SBS21", "SBS26")
pole_sigs_v2  <- c("SBS10", "SBS14")

mmr_sigs_v32  <- c("SBS6", "SBS15", "SBS20", "SBS21", "SBS26", "SBS44")
pole_sigs_v32 <- c("SBS10a", "SBS10b", "SBS14", "SBS28")

get_pathway_scores <- function(df, cosmic_version) {
  if (grepl("v3.2", cosmic_version)) {
    mmr_sigs  <- mmr_sigs_v32
    pole_sigs <- pole_sigs_v32
  } else {
    mmr_sigs  <- mmr_sigs_v2
    pole_sigs <- pole_sigs_v2
  }
  
  df %>%
    mutate(
      pathway = case_when(
        signature %in% mmr_sigs  ~ "MMR",
        signature %in% pole_sigs ~ "POLE",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(pathway)) %>%
    group_by(tool, cosmic_version, dataset, sample, pathway) %>%
    summarise(pathway_score = sum(fraction, na.rm = TRUE), .groups = "drop")
}

pathway_scores <- contrib_df %>%
  group_split(cosmic_version) %>%
  map_dfr(~ get_pathway_scores(.x, unique(.x$cosmic_version)))  # split per CV

# Median per dataset/tool/pathway
pathway_summary <- pathway_scores %>%
  group_by(tool, cosmic_version, dataset, pathway) %>%
  summarise(median_score = median(pathway_score, na.rm = TRUE),
            .groups = "drop")

write_tsv(pathway_scores,
          file.path(out_root, "DeltaMMR_POLE_perSample_referenceOnly.tsv"))
write_tsv(pathway_summary,
          file.path(out_root, "DeltaMMR_POLE_summary_referenceOnly.tsv"))

# Simple plot: MMR / POLE per dataset
p_mmr_pole <- ggplot(pathway_summary,
                     aes(x = dataset, y = median_score,
                         fill = pathway)) +
  geom_col(position = position_dodge(width = 0.7)) +
  facet_grid(tool ~ cosmic_version) +
  theme_bw() +
  labs(x = "Dataset (preclinical)",
       y = "Median pathway contribution (fraction)",
       title = "MMR / POLE signature burden (reference-only)")

ggsave(file.path(out_root, "Figure_DeltaMMR_POLE_referenceOnly.png"),
       p_mmr_pole, width = 9, height = 6, dpi = 300)

msg("  -> ΔMMR / ΔPOLE summaries and plot saved.")

## ------------------------------------------------------------
## 4. WGS vs WES vs TSO (coding vs whole genome approximation)
## ------------------------------------------------------------

msg("Computing WGS vs WES vs TSO comparison (coding vs whole genome)...")

# Here we simply compare datasets directly:
#  - WES + TSO-500 ~ coding
#  - WGS ~ whole genome

# Use total "oncogenic" pathways as examples (MMR + POLE)
pathway_totals <- pathway_scores %>%
  group_by(tool, cosmic_version, dataset, sample) %>%
  summarise(
    MMR  = sum(pathway_score[pathway == "MMR"],  na.rm = TRUE),
    POLE = sum(pathway_score[pathway == "POLE"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = c(MMR, POLE),
               names_to = "pathway",
               values_to = "pathway_score")

write_tsv(pathway_totals,
          file.path(out_root, "WGS_WES_TSO_pathwayTotals_referenceOnly.tsv"))

# Plot: boxplot of pathway_score by dataset
p_coding_vs_whole <- ggplot(pathway_totals,
                            aes(x = dataset, y = pathway_score,
                                fill = dataset)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  facet_grid(tool + pathway ~ cosmic_version) +
  theme_bw() +
  labs(x = "Dataset (WES ~ coding, TSO-500 ~ coding, WGS ~ whole genome)",
       y = "Pathway score (fraction)",
       title = "Coding vs whole genome approximation via WES/TSO-500 vs WGS")

ggsave(file.path(out_root, "Figure_Coding_vs_WholeGenome_Pathways.png"),
       p_coding_vs_whole, width = 10, height = 7, dpi = 300)

msg("  -> WGS vs WES vs TSO coding/whole-genome comparison done.")

## ------------------------------------------------------------
## 5. Flat signature inflation test (WGS vs WES)
## ------------------------------------------------------------

msg("Running flat signature inflation test (WGS vs WES)...")

# Define "flat" signatures – you can expand this set if needed
flat_sigs_v2  <- c("SBS5", "SBS40")
flat_sigs_v32 <- c("SBS5", "SBS40")  # in v3.2 they are still flat; extend if desired

get_flat_scores <- function(df, cosmic_version) {
  if (grepl("v3.2", cosmic_version)) {
    flats <- flat_sigs_v32
  } else {
    flats <- flat_sigs_v2
  }
  
  df %>%
    filter(signature %in% flats) %>%
    group_by(tool, cosmic_version, dataset, sample) %>%
    summarise(flat_score = sum(fraction, na.rm = TRUE), .groups = "drop")
}

flat_scores <- contrib_df %>%
  group_split(cosmic_version) %>%
  map_dfr(~ get_flat_scores(.x, unique(.x$cosmic_version)))

write_tsv(flat_scores,
          file.path(out_root, "FlatSignature_scores_WES_WGS_TSO.tsv"))

# Focus WES vs WGS for inflation test
flat_WES_WGS <- flat_scores %>%
  filter(dataset %in% c("WES", "WGS"))

# Plot
p_flat <- ggplot(flat_WES_WGS,
                 aes(x = dataset, y = flat_score, fill = dataset)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  facet_grid(tool ~ cosmic_version) +
  theme_bw() +
  labs(x = "Dataset",
       y = "Flat signature burden (fraction)",
       title = "Flat signature inflation: WGS vs WES")

ggsave(file.path(out_root, "Figure_FlatSignatureInflation_WGS_vs_WES.png"),
       p_flat, width = 8, height = 6, dpi = 300)

# Simple paired-ish comparison summary (medians)
flat_summary <- flat_WES_WGS %>%
  group_by(tool, cosmic_version, dataset) %>%
  summarise(median_flat = median(flat_score, na.rm = TRUE),
            .groups = "drop")

write_tsv(flat_summary,
          file.path(out_root, "FlatSignatureInflation_summary.tsv"))

msg("  -> Flat signature inflation analysis completed.")

## ------------------------------------------------------------
## Done
## ------------------------------------------------------------

msg("All analyses finished. Outputs in:", out_root)

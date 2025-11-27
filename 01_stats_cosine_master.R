#!/usr/bin/env Rscript

library(tidyverse)
library(ARTool)

## -------------------- Paths --------------------
cos_file   <- "Output/cosine_sim_results/CosineResults_AllTools_AllCOSMIC.tsv"
plotdir    <- "Output/cosine_plots"
resdir     <- "Output/cosine_stats"
mut_outdir <- "Output/mutation_vs_cosine"

dir.create(plotdir,    showWarnings = FALSE, recursive = TRUE)
dir.create(resdir,     showWarnings = FALSE, recursive = TRUE)
dir.create(mut_outdir, showWarnings = FALSE, recursive = TRUE)

## -------------------- Load cosine results --------------------
cosine_all <- read.table(
  cos_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE        # keep "TSO-500" etc.
)

# Expect at least: Sample, Data_Type, Tool, COSMIC_Version, Cosine_Similarity


## ============================================================
## 1) GLOBAL ANOVA + BOXPLOTS
## ============================================================

# --- Boxplot overview ---
p_box <- ggplot(cosine_all,
                aes(x = Data_Type, y = Cosine_Similarity, fill = Tool)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  facet_wrap(~ COSMIC_Version) +
  theme_bw() +
  labs(
    title = "Cosine similarity across workflows, tools, and COSMIC versions",
    x = "Sequencing workflow",
    y = "Cosine similarity"
  )

ggsave(
  file.path(plotdir, "Figure_Cosine_Boxplot.png"),
  p_box,
  width = 8, height = 5, dpi = 300
)

# --- Non-parametric factorial ANOVA (ARTool) ---
cosine_clean <- cosine_all %>%
  drop_na(Cosine_Similarity, Data_Type, Tool, COSMIC_Version) %>%
  mutate(
    Data_Type      = factor(Data_Type),
    Tool           = factor(Tool),
    COSMIC_Version = factor(COSMIC_Version)
  )

art_res   <- art(Cosine_Similarity ~ Data_Type * Tool * COSMIC_Version,
                 data = cosine_clean)
anova_res <- anova(art_res)

# Save ANOVA table
capture.output(
  anova_res,
  file = file.path(resdir, "ART_ANOVA_Cosine.txt")
)


## ============================================================
## 2) PAIRED WILCOXON (WGS vs WES vs TSO-500) + Δ-COS PLOTS
## ============================================================

# Wide format: each row = Sample + Tool + COSMIC_Version with 3 cols: WGS, WES, TSO-500
wide <- cosine_all %>%
  select(Sample, Data_Type, Tool, COSMIC_Version, Cosine_Similarity) %>%
  tidyr::pivot_wider(
    names_from  = Data_Type,
    values_from = Cosine_Similarity
  )

# Only rows with all three workflows present
triplets <- wide %>%
  filter(!is.na(WGS), !is.na(WES), !is.na(`TSO-500`))

# Paired Wilcoxon per Tool × COSMIC_Version
pairwise_p <- triplets %>%
  group_by(Tool, COSMIC_Version) %>%
  summarise(
    p_WGS_vs_WES = wilcox.test(WGS, WES,       paired = TRUE)$p.value,
    p_WES_vs_TSO = wilcox.test(WES, `TSO-500`, paired = TRUE)$p.value,
    p_WGS_vs_TSO = wilcox.test(WGS, `TSO-500`, paired = TRUE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p_WGS_vs_WES_adj = p.adjust(p_WGS_vs_WES, method = "BH"),
    p_WES_vs_TSO_adj = p.adjust(p_WES_vs_TSO, method = "BH"),
    p_WGS_vs_TSO_adj = p.adjust(p_WGS_vs_TSO, method = "BH")
  )

write.table(
  pairwise_p,
  file.path(resdir, "Pairwise_PairedWilcoxon_WGS_WES_TSO.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Δ-cosine per sample
diff_long <- triplets %>%
  mutate(
    diff_WGS_WES = WGS - WES,
    diff_WES_TSO = WES - `TSO-500`,
    diff_WGS_TSO = WGS - `TSO-500`
  ) %>%
  select(Sample, Tool, COSMIC_Version,
         diff_WGS_WES, diff_WES_TSO, diff_WGS_TSO) %>%
  pivot_longer(
    cols      = starts_with("diff_"),
    names_to  = "Comparison",
    values_to = "Delta_Cosine"
  ) %>%
  mutate(
    Comparison = factor(
      Comparison,
      levels = c("diff_WGS_WES", "diff_WES_TSO", "diff_WGS_TSO"),
      labels = c("WGS - WES", "WES - TSO-500", "WGS - TSO-500")
    )
  )

p_diff <- ggplot(diff_long,
                 aes(x = Comparison, y = Delta_Cosine)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(outlier.shape = NA, fill = "grey80") +
  facet_grid(COSMIC_Version ~ Tool) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Pairwise differences in cosine similarity between workflows",
    x = "Comparison",
    y = "Δ Cosine (first - second)"
  )

ggsave(
  file.path(plotdir, "Figure_Cosine_Differences_WilcoxonView.png"),
  p_diff,
  width = 11, height = 7, dpi = 300
)


## ============================================================
## 3) MUTATION COUNT vs COSINE + SPEARMAN CORR
## ============================================================

# --- Load mutation counts (from raw catalogues) ---
load_counts <- function(file) {
  m <- read.table(file, header = TRUE, sep = "\t",
                  row.names = 1, check.names = FALSE)
  tibble(
    Sample    = colnames(m),
    Mutations = colSums(m)
  )
}

counts_WGS <- load_counts("Input files/Preclinical_Dataset_WGS.txt") %>%
  mutate(Data_Type = "WGS")
counts_WES <- load_counts("Input files/Preclinical_Dataset_WES.txt") %>%
  mutate(Data_Type = "WES")
counts_TSO <- load_counts("Input files/Preclinical_Dataset_TSO-500.txt") %>%
  mutate(Data_Type = "TSO-500")

counts_all <- bind_rows(counts_WGS, counts_WES, counts_TSO)

# --- Merge cosine + mutation counts ---
cos_mut <- cosine_all %>%
  left_join(counts_all, by = c("Sample", "Data_Type")) %>%
  filter(!is.na(Mutations), !is.na(Cosine_Similarity))

write.table(
  cos_mut,
  file.path(mut_outdir, "MutCount_vs_Cosine.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# --- Spearman correlation (per Data_Type × Tool × COSMIC_Version) ---
spearman_tbl <- cos_mut %>%
  group_by(Data_Type, Tool, COSMIC_Version) %>%
  summarise(
    n       = n(),
    rho     = suppressWarnings(
      cor(Mutations, Cosine_Similarity, method = "spearman")
    ),
    p_value = suppressWarnings(
      cor.test(Mutations, Cosine_Similarity,
               method = "spearman")$p.value
    ),
    .groups = "drop"
  )

write.table(
  spearman_tbl,
  file.path(mut_outdir, "Spearman_MutCount_vs_Cosine.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# --- Plot mutation count vs cosine similarity (log10 x, dashed 0.9) ---
p_mut <- ggplot(cos_mut,
                aes(x = Mutations,
                    y = Cosine_Similarity,
                    colour = Data_Type)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_hline(yintercept = 0.9, linetype = "dashed") +
  geom_smooth(se = TRUE, method = "loess", linewidth = 0.8) +
  scale_x_log10() +
  theme_bw() +
  labs(
    title = "Mutation count vs cosine similarity",
    x = "Mutation count (log10)",
    y = "Cosine similarity"
  )

ggsave(
  file.path(mut_outdir, "MutationCount_vs_Cosine.png"),
  p_mut,
  width = 7, height = 5, dpi = 300
)


message("✅ 01_stats_cosine_master.R finished.")

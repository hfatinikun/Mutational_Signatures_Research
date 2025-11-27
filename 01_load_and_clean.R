#!/usr/bin/env Rscript

library(tidyverse)

# === Load ===
cosine_all <- read.table(
  "cosine_sim_results/CosineResults_AllTools_AllCOSMIC.tsv",
  header = TRUE, sep = "\t", stringsAsFactors = FALSE
)

# === Clean factors ===
cosine_all <- cosine_all %>%
  mutate(
    Data_Type      = factor(Data_Type, levels = c("WGS", "WES", "TSO-500")),
    Tool           = factor(Tool),
    COSMIC_Version = factor(COSMIC_Version)
  )

# === Optional: add mutation counts ===
source("utils.R")  # for load_crc_data()

get_counts <- function(dataset_type) {
  m <- load_crc_data(dataset_type, input_dir, base_dir)
  tibble(Sample = colnames(m),
         Mutation_Count = colSums(m),
         Data_Type = dataset_type)
}

mut_counts <- bind_rows(
  get_counts("WGS"),
  get_counts("WES"),
  get_counts("TSO-500")
)

cosine_all <- left_join(cosine_all, mut_counts,
                        by = c("Sample", "Data_Type"))

# === Save cleaned data for next scripts ===
write.table(cosine_all,
            "cosine_sim_results/cosine_all_cleaned.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

message("✅ Cleaned cosine data saved.")



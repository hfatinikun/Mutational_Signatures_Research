library(tidyverse)

# Define the folder where the files are located
folder <- "Cosine_Sim_Results"

# List all matching files within that folder
files <- list.files(
  path = folder,
  pattern = "^CosineResults_.*_ALL.tsv$",
  full.names = TRUE
)

# Read and combine all matching files
cosine_all <- files |>
  map(~ read.table(.x, header = TRUE, sep = "\t", stringsAsFactors = FALSE)) |>
  bind_rows()

# Convert columns to factors
cosine_all$Data_Type      <- factor(cosine_all$Data_Type, levels = c("WGS", "WES", "TSO-500"))
cosine_all$Tool           <- factor(cosine_all$Tool)
cosine_all$COSMIC_Version <- factor(cosine_all$COSMIC_Version)

# Write combined results to file (in current directory)
write.table(
  cosine_all,
  "CosineResults_AllTools_AllCOSMIC.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

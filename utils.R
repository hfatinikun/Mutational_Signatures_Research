## ---------- utils.R ----------
## Small helper library for mutational-signature workflows
## Dependencies used by main Rscripts
#cat("\014")

library(MutationalPatterns)# fitting + cosine
library(deconstructSigs)
library(signature.tools.lib)
library(siglasso)
library(ggplot2)             # plotting
library(ggh4x)               # panel size helper
library(grid)                # unit()
library(dplyr)

#Name of files and directory for inputs
base_dir   <- getwd()
cosmic_dir <- "Cosmic files"
input_dir  <- "Input files"

#Function to load COSMIC files
load_cosmic_signatures <- function(cosmic_file, cosmic_dir, base_dir) {
  # Create full path using file.path() for safety (handles spaces and platform differences)
  full_path <- file.path(base_dir, cosmic_dir, cosmic_file)
  # Read the file
  cosmic_data <- read.table(full_path, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
  # Convert to numeric matrix
  cosmic_data <- as.matrix(cosmic_data)
  return(cosmic_data)
}


#Function to load CRC files
load_crc_data <- function(dataset_type = c("WES", "WGS", "TSO-500"),
                          input_dir, base_dir) {
  # Match the dataset type to the correct filename
  dataset_type <- match.arg(dataset_type)  # Ensures valid input
  file_map <- list(
    WES = "Preclinical_Dataset_WES.txt",
    WGS = "Preclinical_Dataset_WGS.txt",
    `TSO-500`  = "Preclinical_Dataset_TSO-500.txt"
  )
  # Build the full path to the file
  full_path <- file.path(base_dir, input_dir, file_map[[dataset_type]])
  # Read the file
  crc_data <- read.table(full_path, header = TRUE, row.names = 1)
  return(crc_data)
}


#Function to fit CRC data with MutationalPatterns
fit_crc_with_MutationalPatterns <- function(crc_data, cosmic_data, dataset_type = c("WES", "WGS", "TSO-500")) {
  dataset_type <- match.arg(dataset_type)
  print(paste("[MutationalPatterns] Fitting", dataset_type))
  
  fitting_result <- fit_to_signatures(crc_data, cosmic_data)
  return(fitting_result)
}


#Function to normalize the contributions
normalize_signature_contributions <- function(fitting_result) {
  contribution <- fitting_result$contribution
  normalized <- contribution / colSums(contribution)[col(contribution)]
  return(as.data.frame(normalized))
}


#Function for Cosine similarity Analysis
calculate_cosine_similarity <- function(fitting_result_WES, fitting_result_TSO_500, fitting_result_WGS,
                                        CRC_WES, CRC_TSO_500, CRC_WGS) {
  
  # Calculate cosine similarities for each dataset
  cosine_sim_CRC_WES_cosmic_MP <- diag(cos_sim_matrix(fitting_result_WES$reconstructed, CRC_WES))
  cosine_sim_CRC_TSO500_cosmic_MP <- diag(cos_sim_matrix(fitting_result_TSO_500$reconstructed, CRC_TSO_500))
  cosine_sim_CRC_WGS_cosmic_MP    <- diag(cos_sim_matrix(fitting_result_WGS$reconstructed, CRC_WGS))
  
  # Create data frames for each dataset type
  cosine_data_type_WES <- data.frame(cosine_sim_CRC_WES_cosmic_MP, "WES")
  colnames(cosine_data_type_WES) <- c("Cosine_Similarity", "Data_Type")
  
  cosine_data_type_TSO500 <- data.frame(cosine_sim_CRC_TSO500_cosmic_MP, "TSO-500")
  colnames(cosine_data_type_TSO500) <- c("Cosine_Similarity", "Data_Type")
  
  cosine_data_type_WGS <- data.frame(cosine_sim_CRC_WGS_cosmic_MP, "WGS")
  colnames(cosine_data_type_WGS) <- c("Cosine_Similarity", "Data_Type")
  
  # Combine all the data frames into one
  cosine_data_type <- rbind(cosine_data_type_WES, cosine_data_type_TSO500, cosine_data_type_WGS)
  colnames(cosine_data_type) <- c("Cosine_Similarity", "Data_Type")
  
  return(cosine_data_type)
}

# Extended cosine similarity function for merging results
# Improved version with Sample column
calculate_cosine_similarity_extended <- function(fitting_result_WES, fitting_result_TSO_500, fitting_result_WGS,
                                                 CRC_WES, CRC_TSO_500, CRC_WGS, tool_name, cosmic_version) {
  # --- Calculate cosine similarities and keep sample names
  cosine_WES <- diag(cos_sim_matrix(fitting_result_WES$reconstructed, CRC_WES))
  cosine_TSO <- diag(cos_sim_matrix(fitting_result_TSO_500$reconstructed, CRC_TSO_500))
  cosine_WGS <- diag(cos_sim_matrix(fitting_result_WGS$reconstructed, CRC_WGS))
  
  # --- Create tidy data frames for each dataset
  df_WES <- data.frame(
    Sample = names(cosine_WES),
    Cosine_Similarity = as.numeric(cosine_WES),
    Data_Type = "WES",
    Tool = tool_name,
    COSMIC_Version = cosmic_version
  )
  
  df_TSO <- data.frame(
    Sample = names(cosine_TSO),
    Cosine_Similarity = as.numeric(cosine_TSO),
    Data_Type = "TSO-500",
    Tool = tool_name,
    COSMIC_Version = cosmic_version
  )
  
  df_WGS <- data.frame(
    Sample = names(cosine_WGS),
    Cosine_Similarity = as.numeric(cosine_WGS),
    Data_Type = "WGS",
    Tool = tool_name,
    COSMIC_Version = cosmic_version
  )
  
  # --- Combine all datasets into one
  cosine_data <- rbind(df_WGS, df_WES, df_TSO)
  
  return(cosine_data)
}



#Function to Create the Boxplot for Cosine Similarity
plot_cosine_similarity <- function(cosine_data, title = "Cosine Similarity with Cosmic v3.2_SBS_GRCh38(SLT)") {
  
  # Create the ggplot for the boxplot
  p <- ggplot(cosine_data,
              aes(x = factor(Data_Type, levels = c("WGS", "WES", "TSO-500")),
                  y = Cosine_Similarity)) +
    geom_boxplot(color = "black", fill = "grey") +
    coord_flip() +
    labs(x = "NGS Workflow", y = "Cosine Similarity") +
    ggtitle(title) +
    scale_y_continuous(breaks = seq(0.9, 1, 0.025),
                       limits = c(0.89, 1),
                       expand = c(0, 0)) +
    geom_hline(yintercept = 0.9, linetype = "dashed") +
    force_panelsizes(rows = unit(6, "cm"), cols = unit(8, "cm")) +
    theme_bw() +
    theme(axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          axis.text     = element_text(size = 10),
          panel.border  = element_rect(color = "black", linewidth = 0.4))
  
  return(p)
}


#--------------DeconstructSigs-------------------------

#Function to fit CRC data with deconstructSigs
# Fit CRC data with deconstructSigs and return contributions + reconstructed profiles
fit_crc_with_deconstructSigs <- function(crc_data, cosmic_data, dataset_type = c("WES", "WGS", "TSO-500")) {
  dataset_type <- match.arg(dataset_type)
  print(paste("[deconstructSigs] Fitting", dataset_type))
  
  crc_t    <- as.data.frame(t(crc_data))      # samples × contexts
  cosmic_t <- as.data.frame(t(cosmic_data))   # signatures × contexts
  common   <- intersect(colnames(crc_t), colnames(cosmic_t))
  crc_t    <- crc_t[, common]
  cosmic_t <- cosmic_t[, common]
  
  res_list <- lapply(rownames(crc_t), function(s) {
    cat("Processing sample", s, "\n")
    out <- deconstructSigs::whichSignatures(
      tumor.ref        = crc_t,
      sample.id        = s,
      signatures.ref   = cosmic_t,
      contexts.needed  = TRUE,
      signature.cutoff = 0
    )
    w <- out$weights[1, , drop = FALSE]
    rownames(w) <- s
    w
  })
  
  weights <- do.call(rbind, res_list)         # samples × signatures
  contrib <- t(weights)                       # signatures × samples
  recon   <- t(cosmic_t) %*% contrib          # contexts × samples
  
  list(contribution = contrib,
       reconstructed = recon,
       method = "deconstructSigs",
       dataset_type = dataset_type)
}

#----------------------------signature.lib.tools-----------------------------------
#Function to fit CRC data with signature.lib.tools
fit_crc_with_signature_tools_lib <- function(crc_data, cosmic_data, dataset_type = c("WES", "WGS", "TSO-500")) {
  dataset_type <- match.arg(dataset_type)
  print(paste("[signature.tools.lib] Fitting", dataset_type))
  
  # crc_data: contexts × samples
  # cosmic_data: contexts × signatures
  res <- signature.tools.lib::Fit(
    catalogues          = as.matrix(crc_data),
    signatures          = as.matrix(cosmic_data),
    exposureFilterType  = "giniScaledThreshold",
    useBootstrap        = TRUE,
    nboot               = 100,
    nparallel           = 1
  )
  
  # Exposures -> contributions (signatures × samples), drop "unassigned"
  w <- t(res$exposures)  # signatures × samples
  if ("unassigned" %in% rownames(w)) {
    w <- w[setdiff(rownames(w), "unassigned"), , drop = FALSE]
  }
  
  # Ensure the signature order matches between cosmic_data (columns) and w (rows)
  common_sigs <- intersect(colnames(cosmic_data), rownames(w))
  if (length(common_sigs) == 0) {
    stop("No overlapping signatures between cosmic_data and exposures.")
  }
  cosmic_use <- as.matrix(cosmic_data[, common_sigs, drop = FALSE])      # contexts × signatures
  contribution <- as.matrix(w[common_sigs, , drop = FALSE])              # signatures × samples
  
  # Reconstruct catalogues (contexts × samples)
  reconstructed <- cosmic_use %*% contribution
  
  list(
    contribution  = contribution,      # signatures × samples
    reconstructed = reconstructed,     # contexts × samples
    method        = "signature.tools.lib",
    dataset_type  = dataset_type
  )
}


#----------------sigLASSO----------------
fit_crc_with_sigLASSO <- function(crc_data, cosmic_data, dataset_type = c("WES", "WGS", "TSO-500")) {
  dataset_type <- match.arg(dataset_type)
  message(paste("[sigLASSO] Fitting", dataset_type))
  
  # Ensure matrices
  spectrum  <- as.matrix(crc_data)    # contexts × samples
  signature <- as.matrix(cosmic_data) # contexts × signatures
  
  # 1) Align mutation contexts (rows)
  common_ctx <- intersect(rownames(spectrum), rownames(signature))
  if (length(common_ctx) == 0L) stop("No overlapping mutation contexts between CRC and COSMIC.")
  spectrum  <- spectrum [common_ctx, , drop = FALSE]
  signature <- signature[common_ctx, , drop = FALSE]
  
  # 2) Run sigLASSO
  # (set.seed if you want reproducibility)
  res <- siglasso::siglasso(spectrum, signature = signature)  # returns signatures × samples
  
  # 3) Coerce to numeric matrix (no transpose!)
  contribution <- as.matrix(res)                # signatures × samples
  
  # 4) Align signatures by name (defensive)
  common_sigs <- intersect(colnames(signature), rownames(contribution))
  if (length(common_sigs) == 0L) stop("No overlapping signature names between result and COSMIC.")
  signature    <- signature[, common_sigs, drop = FALSE]
  contribution <- contribution[common_sigs, , drop = FALSE]
  
  # 5) Reconstruct catalogues (contexts × samples)
  reconstructed <- signature %*% contribution
  
  list(
    contribution  = contribution,   # signatures × samples
    reconstructed = reconstructed,  # contexts × samples
    method        = "sigLASSO",
    dataset_type  = dataset_type
  )
}
  
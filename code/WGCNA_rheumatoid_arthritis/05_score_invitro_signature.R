# # ==================================================================
# SCRIPT: Calculate Patient-Level Gene Set Scores and Deltas
#
# Description:
# - Loads the 90 gene sets from "IM049_top100_up.csv".
# - Calculates the average LogCPM expression score for each gene set
#   in every paired TCZ patient sample (W0 and W16).
# - Calculates the Delta score (W16 - W0) for each patient/gene set pair.
#
# Outputs:
# - IM049_Top100_GeneSet_Deltas_TCZ.csv
# ==================================================================

# --- Step 1: Load Libraries and Prepare IM049 Gene Sets ---
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(readr)

# Load the gene sets file created by the DEA script
IM049_DEG_FILE <- "IM049_top100_up.csv"
im049_data <- read.csv(IM049_DEG_FILE, check.names = FALSE)

# Create the final list format: {contrast_name: [genes]}
im049_gene_sets <- im049_data %>%
  dplyr::select(contrast, gene_name) %>%
  distinct() %>%
  split(.$contrast) %>%
  lapply(function(x) x$gene_name)

num_im049_sets <- length(im049_gene_sets)
cat(paste0("### Loaded ", num_im049_sets, " IM049 gene sets. ###\n"))


# --- Step 2: Define and Run Scoring for TCZ (TOC) ---

drug_filter <- "TOC"
metadata_file <- "R4RA_metadata.csv"
expr_file <- "TCZ_data_proteincoding.csv"
output_file <- "IM049_Top100_GeneSet_Deltas_TCZ.csv"

cat(paste0("\n### Calculating Average logCPM for ", drug_filter, " cohort... ###\n"))

# Load and filter metadata (W0, W16 for TOC)
metadata_raw <- read.csv(metadata_file)
pheno_subset <- metadata_raw %>% filter(Drug == drug_filter, Timepoint %in% c(0, 16))

# Load expression data
datExpr <- read.csv(expr_file) %>%
  tibble::column_to_rownames("gene_name") %>%
  as.matrix()
colnames(datExpr) <- str_replace(colnames(datExpr), "^X", "")

datExpr_aligned <- datExpr[, colnames(datExpr) %in% pheno_subset$Header, drop = FALSE]

if (ncol(datExpr_aligned) == 0) {
  stop(paste("No expression data found for", drug_filter, "W0/W16 samples. Aborting."))
}

score_list <- list()
paired_patient_ids <- pheno_subset %>% 
  group_by(Patient_ID) %>% 
  filter(n() == 2) %>% 
  pull(Patient_ID) %>% 
  unique()

pheno_subset_paired <- pheno_subset %>% filter(Patient_ID %in% paired_patient_ids)

for (set_name in names(im049_gene_sets)) {
  gene_set <- im049_gene_sets[[set_name]]
  
  genes_in_data <- intersect(gene_set, rownames(datExpr_aligned))
  
  if (length(genes_in_data) >= 5) {
    subset_expr <- datExpr_aligned[genes_in_data, , drop = FALSE]
    
    # Calculate the average logCPM score for all paired samples
    avg_logcpm_scores <- colMeans(subset_expr, na.rm = TRUE)
    
    score_list[[set_name]] <- data.frame(
      Header = names(avg_logcpm_scores),
      Gene_Set = set_name,
      Avg_logCPM = avg_logcpm_scores
    )
  }
}

# --- Step 3: Calculate Paired Delta and Save ---

tcz_scores <- bind_rows(score_list)

final_summary_df_tcz <- tcz_scores %>%
  # Join with paired metadata subset
  left_join(pheno_subset_paired %>% dplyr::select(Header, Patient_ID, Timepoint, Drug), by = "Header") %>%
  dplyr::select(Patient_ID, Drug, Gene_Set, Timepoint, Avg_logCPM) %>%
  tidyr::drop_na() %>%
  # Pivot wider for paired delta calculation
  tidyr::pivot_wider(
    names_from = Timepoint,
    values_from = Avg_logCPM,
    names_prefix = "W"
  ) %>%
  # Calculate delta score
  mutate(Delta_logCPM = W16 - W0) %>%
  # Ensure only fully paired samples are kept (already mostly done, but confirms structure)
  tidyr::drop_na(W0, W16)

write.csv(final_summary_df_tcz, output_file, row.names = FALSE)
cat(paste0("\nTCZ patient scoring complete. Results saved to '", output_file, "'\n"))

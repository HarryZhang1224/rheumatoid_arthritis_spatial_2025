# ==================================================================
# R SCRIPT: All-vs-All Differential Expression Analysis (90 UP Signatures)
#
# Description:
# - Generates 90 unique directional gene signatures by running all 90 non-identical
#   pairwise contrasts (A vs B and B vs A).
# - Extracts the top 100 UP-regulated genes (LFC > 0) from each contrast,
#   filtered by an unadjusted P-value cutoff of 0.05.
#
# Outputs:
# - IM049_top100_up.csv
# ==================================================================

# --- Step 1: Load Libraries ---
library(limma)
library(dplyr)
library(tibble)
library(stringr)

# --- Step 2: Define Files ---
EXPRESSION_FILE <- "IM049_data_proteincoding.csv"
METADATA_FILE <- "IM049_metadata.csv"
OUTPUT_FILE <- "IM049_top100_up.csv"

GENE_COL <- "gene_name"
SAMPLE_ID_COL <- "SampleID"
P_VALUE_CUTOFF <- 0.05

cat("### Loading expression and metadata files... ###\n")

# --- Step 3: Load Data ---
expr_data <- read.csv(EXPRESSION_FILE, check.names = FALSE)
metadata <- read.csv(METADATA_FILE, check.names = FALSE)

# --- Step 4: Fix sample names and ensure common samples ---
expr_samples <- colnames(expr_data)[colnames(expr_data) != GENE_COL]
fixed_expr_samples <- str_replace_all(expr_samples, "^X", "")
colnames(expr_data)[colnames(expr_data) != GENE_COL] <- fixed_expr_samples

# Align metadata and expression data
common_samples <- intersect(fixed_expr_samples, metadata[[SAMPLE_ID_COL]])
expr_data <- expr_data[, c(GENE_COL, common_samples)]
metadata <- metadata[metadata[[SAMPLE_ID_COL]] %in% common_samples, ]

# --- Step 5: Prepare expression matrix and conditions ---
expr_matrix <- expr_data %>%
  tibble::column_to_rownames(var = GENE_COL)

# Find unique condition column names (R-safe names created by read.csv)
condition_cols <- setdiff(colnames(metadata), SAMPLE_ID_COL)
metadata_tibble <- as_tibble(metadata)
all_results <- list()

# --- Step 6: Run all-vs-all DEA (90 Contrasts) ---
cat("### Running all 90 directional differential expression contrasts... ###\n")

for (cond1_name in condition_cols) {
  for (cond2_name in condition_cols) {
    if (cond1_name == cond2_name) next # Skip self-comparison (A vs A)
    
    # Build contrast string Cond2 vs Cond1
    contrast_name <- paste(cond2_name, "vs", cond1_name)
    cat("Comparing:", contrast_name, "\n")
    
    # --- Data Subset and Design Matrix ---
    group_cols <- c(cond1_name, cond2_name)
    
    # Select relevant sample columns
    group <- metadata_tibble %>%
      dplyr::select(!!sym(SAMPLE_ID_COL), all_of(group_cols))
    
    # Determine which samples belong to which group
    group_vector <- ifelse(group[[cond1_name]] == 1, cond1_name,
                           ifelse(group[[cond2_name]] == 1, cond2_name, NA))
    
    keep <- !is.na(group_vector)
    group_vector <- group_vector[keep]
    expr_subset <- expr_matrix[, keep]
    
    # Design matrix setup
    design <- model.matrix(~0 + factor(group_vector))
    colnames(design) <- c(cond1_name, cond2_name)
    
    # limma contrast
    contrast <- makeContrasts(contrasts = paste(cond2_name, "-", cond1_name), levels = design)
    fit <- lmFit(expr_subset, design)
    fit2 <- contrasts.fit(fit, contrast)
    fit2 <- eBayes(fit2)
    
    # --- Gene Signature Extraction ---
    
    # Extract all results
    top_table <- topTable(fit2, number = Inf, sort.by = "B") %>%
      tibble::rownames_to_column(var = GENE_COL)
    
    # Filter for UP-regulation and P-value cutoff
    signature_df <- top_table %>%
      filter(P.Value < P_VALUE_CUTOFF) %>%       # Apply P-value cutoff
      filter(logFC > 0) %>%                      # Select only UP-regulated genes
      arrange(desc(logFC)) %>%                   # Sort by LFC magnitude (most positive first)
      head(100) %>%                              # Select top 100
      mutate(contrast = contrast_name, .before = 1) %>% # Renamed to 'contrast' for compatibility
      mutate(direction = "Upregulated")
    
    if (nrow(signature_df) > 0) {
      all_results[[contrast_name]] <- signature_df
    }
  }
}

# --- Step 7: Combine results and save ---
IM049_final_results_df <- bind_rows(all_results)
IM049_final_results_df <- IM049_final_results_df %>%
  dplyr::select(contrast, gene_name, logFC, P.Value, adj.P.Val, direction)

write.csv(IM049_final_results_df, OUTPUT_FILE, row.names = FALSE)
cat("### Differential expression analysis complete. Results saved to", OUTPUT_FILE, "###\n")
cat(paste0("Total rows (genes in all signatures): ", nrow(IM049_final_results_df), "\n"))

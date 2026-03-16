# ==================================================================
# SCRIPT: Prepare Merged Protein-Coding Expression Dataset
#
# Description:
# - STEP A: Filters raw expression files (TCZ, RTX, IM049) to keep only protein-coding genes.
# - STEP B: Merges the protein-coding filtered TCZ and RTX files.
#
# Inputs:
# - Raw expression CSV files (specified in 'RAW_EXPRESSION_FILES' list below).
#   Each file should have genes in rows, samples in columns, and a column
#   (assumed 'gene_name' or first column) containing gene symbols.
#
# Intermediates (created by Step A, used by Step B):
# - *_proteincoding.csv files (e.g., "TCZ_data_proteincoding.csv", "RTX_data_proteincoding.csv")
#
# Final Output:
# - Merged_data_proteincoding.csv: Combined logCPM matrix (Genes x Samples)
#                                 containing only common protein-coding genes
#                                 and all TCZ/RTX samples.
# ==================================================================

# --- Load Libraries ---
# Install BiocManager and biomaRt if needed
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# if (!requireNamespace("biomaRt", quietly = TRUE)) BiocManager::install("biomaRt")
# if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
# if (!requireNamespace("tibble", quietly = TRUE)) install.packages("tibble")

library(biomaRt)
library(dplyr)
library(tibble)

# ==================================================================
# STEP A: Filter Raw Files to Protein-Coding Genes Only
# ==================================================================

# --- Define Raw Input Files ---
# !!! USER ACTION: Modify this list to include your ACTUAL raw input files !!!
RAW_EXPRESSION_FILES <- c("TCZ_data.csv", "RTX_data.csv", "IM049_data.csv") # Add others as needed

# Define filenames for the protein-coding filtered outputs (used later)
TCZ_PC_FILE <- sub("\\.csv$", "_proteincoding.csv", "TCZ_data.csv")
RTX_PC_FILE <- sub("\\.csv$", "_proteincoding.csv", "RTX_data.csv")

# --- Get Protein-Coding Gene List from Ensembl ---
cat("### STEP A: Filtering Raw Files to Protein-Coding Genes ###\n")
cat(" -> Retrieving protein-coding gene list from Ensembl...\n")
protein_coding_genes <- c() # Initialize
tryCatch({
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  gene_info <- getBM(attributes = c('hgnc_symbol', 'gene_biotype'), mart = ensembl)
  
  protein_coding_genes <- unique(gene_info$hgnc_symbol[gene_info$gene_biotype == "protein_coding"])
  protein_coding_genes <- protein_coding_genes[!is.na(protein_coding_genes) & protein_coding_genes != ""]
  
  if (length(protein_coding_genes) == 0) {
    stop("Failed to retrieve a valid list of protein-coding genes from Ensembl.")
  }
  cat(paste0(" -> Retrieved ", length(protein_coding_genes), " unique protein-coding gene symbols.\n"))
  
}, error = function(e){
  stop("Error connecting to or querying Ensembl via biomaRt: ", conditionMessage(e),
       "\nPlease check your internet connection and Ensembl status.")
})

# --- Read, Filter, and Save Each Raw Data File ---
cat("\n -> Processing and filtering raw input expression files...\n")

filtered_count <- 0
for (file in RAW_EXPRESSION_FILES) {
  if (!file.exists(file)) {
    warning(paste(" -> Input file not found:", file, "- Skipping."))
    next
  }
  
  cat(paste0("  -> Processing: ", file, "\n"))
  original_df <- read.csv(file, check.names = FALSE)
  
  gene_col_name <- NULL
  if ("gene_name" %in% colnames(original_df)) {
    gene_col_name <- "gene_name"
  } else if (is.character(original_df[[1]])) {
    gene_col_name <- colnames(original_df)[1]
    cat(paste0("     -> Using first column '", gene_col_name, "' as gene names.\n"))
  } else {
    warning(paste0("     -> Cannot identify gene name column in ", file, ". Skipping filtering."))
    next
  }
  
  original_rows <- nrow(original_df)
  #cat(paste0("     -> Original dimensions: ", original_rows, " rows x ", ncol(original_df), " columns.\n"))
  
  filtered_df <- original_df[original_df[[gene_col_name]] %in% protein_coding_genes, ]
  filtered_rows <- nrow(filtered_df)
  cat(paste0("     -> Filtered dimensions: ", filtered_rows, " rows x ", ncol(filtered_df), " columns.\n"))
  
  if (filtered_rows > 0) {
    new_filename <- sub("\\.csv$", "_proteincoding.csv", file)
    write.csv(filtered_df, file = new_filename, row.names = FALSE)
    cat(paste0("     -> Successfully saved filtered data to: ", new_filename, "\n"))
    filtered_count <- filtered_count + 1
  } else {
    warning(paste0("     -> No protein-coding genes found in ", file, ". No output file saved.\n"))
  }
}

cat(paste0(" -> STEP A complete. Processed ", filtered_count, " files successfully.\n"))


# ==================================================================
# STEP B: Merge Protein-Coding TCZ and RTX Files
# ==================================================================

cat("\n### STEP B: Merging Protein-Coding TCZ and RTX Files ###\n")

MERGED_OUTPUT_FILE <- "Merged_data_proteincoding.csv"

# --- Load Individual Protein-Coding Datasets ---
cat(" -> Loading individual TCZ and RTX protein-coding datasets...\n")

# Check if required input files exist (created by Step A)
if (!file.exists(TCZ_PC_FILE)) stop("ERROR: Protein-coding TCZ file not found: ", TCZ_PC_FILE)
if (!file.exists(RTX_PC_FILE)) stop("ERROR: Protein-coding RTX file not found: ", RTX_PC_FILE)

# Load TCZ data
tcz_pc_data <- read.csv(TCZ_PC_FILE, check.names = FALSE)
if (!"gene_name" %in% colnames(tcz_pc_data)) {
  if (is.character(tcz_pc_data[[1]])) { tcz_pc_data <- tcz_pc_data %>% rename(gene_name = 1) }
  else { stop("ERROR: Cannot identify gene name column in ", TCZ_PC_FILE) }
}
tcz_pc_data <- tcz_pc_data %>% distinct(gene_name, .keep_all = TRUE) %>% column_to_rownames("gene_name")
# cat(paste0(" -> TCZ PC data loaded: ", nrow(tcz_pc_data), " genes x ", ncol(tcz_pc_data), " samples.\n"))

# Load RTX data
rtx_pc_data <- read.csv(RTX_PC_FILE, check.names = FALSE)
if (!"gene_name" %in% colnames(rtx_pc_data)) {
  if (is.character(rtx_pc_data[[1]])) { rtx_pc_data <- rtx_pc_data %>% rename(gene_name = 1) }
  else { stop("ERROR: Cannot identify gene name column in ", RTX_PC_FILE) }
}
rtx_pc_data <- rtx_pc_data %>% distinct(gene_name, .keep_all = TRUE) %>% column_to_rownames("gene_name")
# cat(paste0(" -> RTX PC data loaded: ", nrow(rtx_pc_data), " genes x ", ncol(rtx_pc_data), " samples.\n"))

# --- Identify Common Genes ---
cat(" -> Identifying common protein-coding genes...\n")
common_genes <- intersect(rownames(tcz_pc_data), rownames(rtx_pc_data))

if (length(common_genes) == 0) {
  stop("ERROR: No common protein-coding genes found between the TCZ and RTX datasets.")
}
cat(paste0(" -> Found ", length(common_genes), " common protein-coding genes.\n"))

# --- Filter and Align Datasets ---
cat(" -> Filtering datasets to common genes and ensuring same gene order...\n")
tcz_data_filtered <- tcz_pc_data[common_genes, , drop = FALSE]
rtx_data_filtered <- rtx_pc_data[common_genes, , drop = FALSE] # Already aligned by common_genes index
stopifnot(all(rownames(tcz_data_filtered) == rownames(rtx_data_filtered))) # Sanity check

# --- Merge Datasets ---
cat(" -> Merging TCZ and RTX datasets column-wise...\n")
overlapping_samples <- intersect(colnames(tcz_data_filtered), colnames(rtx_data_filtered))
if (length(overlapping_samples) > 0) {
  warning("WARNING: Overlapping sample names found: ", paste(overlapping_samples, collapse=", "))
}
merged_data <- cbind(tcz_data_filtered, rtx_data_filtered)
cat(paste0(" -> Merged data dimensions: ", nrow(merged_data), " genes x ", ncol(merged_data), " samples.\n"))

# --- Save Merged Dataset ---
cat(" -> Saving merged dataset...\n")
merged_data_to_save <- merged_data %>%
  as.data.frame() %>%
  rownames_to_column("gene_name")
write.csv(merged_data_to_save, MERGED_OUTPUT_FILE, row.names = FALSE)
cat(paste0(" -> Merged data saved successfully to '", MERGED_OUTPUT_FILE, "'.\n"))

cat("\n### Data preparation script complete. ###\n")

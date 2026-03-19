# ==================================================================
# SCRIPT: Run WGCNA on Merged TCZ+RTX Bulk RNA-seq Data
#
# Description:
# - Determines the optimal soft power threshold.
# - Constructs the co-expression network using blockwiseModules.
# - Identifies gene modules.
# - Calculates Module Eigengenes (MEs) and Gene Module Membership (kME).
#
# Inputs:
# - Merged_data_proteincoding.csv: logCPM expression matrix (Genes x Samples)
#
# Outputs:
# - Merged_WGCNA_soft_power_selection.png: Plots for soft power selection.
# - Merged_WGCNA_module_assignments.csv: Gene name, module color, module number.
# - Merged_WGCNA_module_dendrogram.png: Dendrogram with module colors.
# - Merged_WGCNA_gene_kME.csv: Gene name, module color/number, kME value for assigned module.
# - Merged_WGCNA_net_object.RData: R object containing the 'merged_net' results from blockwiseModules.
# ==================================================================

# --- Step 1: Load Libraries ---
# Check and install WGCNA if needed
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# if (!requireNamespace("WGCNA", quietly = TRUE)) BiocManager::install("WGCNA")
# if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
# if (!requireNamespace("tibble", quietly = TRUE)) install.packages("tibble")

library(WGCNA)
library(dplyr)
library(tibble)

# Allow WGCNA to use multiple threads if available
allowWGCNAThreads()

# --- Step 2: Define File Paths ---
EXPR_FILE <- "Merged_data_proteincoding.csv"
OUTPUT_PREFIX <- "Merged_WGCNA"

# --- Step 3: Load and Format Expression Data ---
cat("### Loading and formatting expression data... ###\n")

# Load the protein-coding Merged expression data
# Assumes first column is gene names
Merged_pc_data <- read.csv(EXPR_FILE, row.names = 1)

# Format for WGCNA: Samples as rows, Genes as columns
merged_datExpr <- Merged_pc_data %>%
  as.data.frame() %>%
  t() # Transpose

cat(paste0(" -> Expression matrix dimensions: ", nrow(merged_datExpr), " samples x ", ncol(merged_datExpr), " genes.\n"))

# --- Step 4: Check Data Quality (Genes and Samples) ---
cat("\n### Checking data quality... ###\n")
gsg <- goodSamplesGenes(merged_datExpr, verbose = 3)
if (!gsg$allOK) {
  # Optionally remove problematic genes and samples
  if (sum(!gsg$goodGenes) > 0) {
    cat(paste0(" -> Removing ", sum(!gsg$goodGenes), " genes with too many missing values or zero variance.\n"))
    merged_datExpr <- merged_datExpr[, gsg$goodGenes]
  }
  if (sum(!gsg$goodSamples) > 0) {
    cat(paste0(" -> Removing ", sum(!gsg$goodSamples), " samples with too many missing values.\n"))
    merged_datExpr <- merged_datExpr[gsg$goodSamples, ]
  }
  cat(paste0(" -> New dimensions: ", nrow(merged_datExpr), " samples x ", ncol(merged_datExpr), " genes.\n"))
} else {
  cat(" -> Data quality check passed (all samples and genes OK).\n")
}


# --- Step 5: Select Soft Power Threshold ---
cat("\n### Selecting soft power threshold... ###\n")

# Define a range of powers to test
powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))

# Call the network topology analysis function
# Use tryCatch for robustness if pickSoftThreshold fails
sft <- tryCatch({
  # Note: networkType="unsigned" is needed here to match the network construction type
  # If the original run didn't explicitly set 'networkType' in pickSoftThreshold, 
  # WGCNA defaults may be slightly different. We keep it as "unsigned" as it was 
  # in your first script, which is the standard practice for matching step 6.
  pickSoftThreshold(merged_datExpr, powerVector = powers, verbose = 5, networkType = "unsigned")
}, error = function(e){
  stop("Error in pickSoftThreshold: ", conditionMessage(e),
       "\nCheck input data for NAs or zero variance columns not caught by goodSamplesGenes.")
})


# Find the suggested power (lowest power reaching R^2 > 0.90, or use the one determined previously)
# For this dataset, β = 8 was previously identified.
merged_soft_power <- 8 # Set based on previous analysis or inspection of plots
cat(paste0(" -> Using selected soft power: ", merged_soft_power, "\n"))

# Plot the results
png(paste0(OUTPUT_PREFIX, "_soft_power_selection.png"), width = 900, height = 500)
par(mfrow = c(1, 2))
cex1 <- 0.9 # Scaling factor for text size
# Scale-free topology fit index vs. soft-thresholding power
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n", main = "Scale independence")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red")
abline(h = 0.90, col = "red") # Threshold line (adjust R^2 threshold if needed)
# Mean connectivity vs. soft-thresholding power
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n",
     main = "Mean connectivity")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1, col = "red")
dev.off()
cat(paste0(" -> Soft power selection plots saved to '", OUTPUT_PREFIX, "_soft_power_selection.png'.\n"))

# --- Step 6: Construct Network and Identify Modules ---
cat("\n### Constructing network and identifying modules... ###\n")

# !!! REVERTED BLOCKWISEMODULES CALL TO ORIGINAL PUBLISHED PARAMETERS !!!
# The critical change here is removing 'networkType' and 'maxBlockSize' 
# and changing 'saveTOMFileBase' to match your "used before" snippet.
merged_net <- blockwiseModules(
  merged_datExpr,
  power = merged_soft_power,
  # networkType = "unsigned",  # REMOVED TO MATCH ORIGINAL CALL
  TOMType = "unsigned",      # Use the standard unsigned network type
  minModuleSize = 30,        # Minimum number of genes in a module
  reassignThreshold = 0,     # Standard setting
  mergeCutHeight = 0.25,     # Threshold for merging similar modules
  numericLabels = TRUE,      # Creates numeric labels for modules
  saveTOMs = TRUE,           # Optional: saves the Topological Overlap Matrix
  saveTOMFileBase = "Merged_TOM", # REVERTED TO ORIGINAL FILENAME BASE
  verbose = 3
  # maxBlockSize parameter REMOVED to match original call (relying on WGCNA default)
)

cat(paste0(" -> Found ", length(unique(merged_net$colors)), " modules (including grey module 0).\n"))

# --- Step 7: Save Module Assignments ---
cat("\n### Saving module assignments... ###\n")

# Module number 0 is reserved for unassigned genes (grey module)
module_numbers <- merged_net$colors # Numeric labels are directly in $colors when numericLabels=TRUE
module_colors <- labels2colors(module_numbers) # Get corresponding standard color names

# Create a dataframe with gene names, module numbers, and module colors
merged_module_assignments <- data.frame(
  gene_name = colnames(merged_datExpr), # Genes are columns after transpose
  module_number = module_numbers,
  module_color = module_colors
)

# Save this dataframe to a CSV
write.csv(merged_module_assignments, paste0(OUTPUT_PREFIX, "_module_assignments.csv"), row.names = FALSE)
cat(paste0(" -> Module assignments saved to '", OUTPUT_PREFIX, "_module_assignments.csv'.\n"))

# --- Step 8: Plot and Save the Gene Dendrogram ---
cat("\n### Plotting and saving gene dendrogram... ###\n")

# Check if dendrogram exists (might not if blockwise failed or only one block with issues)
if (length(merged_net$dendrograms) >= 1) {
  png(paste0(OUTPUT_PREFIX, "_module_dendrogram.png"), width = 1200, height = 800)
  # Plot the dendrogram with module colors underneath
  plotDendroAndColors(
    merged_net$dendrograms[[1]], # Assumes only one block was processed
    module_colors[merged_net$blockGenes[[1]]],
    "Module colors",
    dendroLabels = FALSE, # Do not show gene names on dendrogram
    hang = 0.03,
    addGuide = TRUE,
    guideHang = 0.05,
    main = "Gene Dendrogram and Module Colors (Merged Data)"
  )
  dev.off()
  cat(paste0(" -> Dendrogram plot saved to '", OUTPUT_PREFIX, "_module_dendrogram.png'.\n"))
} else {
  cat(" -> Skipping dendrogram plot: No dendrogram found in WGCNA output.\n")
}


# --- Step 9: Calculate and Save Gene Module Membership (kME) ---
cat("\n### Calculating and saving Gene Module Membership (kME) for all genes including grey... ###\n")

# 1. Module Eigengenes (MEs)
MEs <- merged_net$MEs
if (is.null(rownames(MEs))) rownames(MEs) <- rownames(merged_datExpr)

# 2. Compute kME: correlation of each gene with all module eigengenes
gene_kME <- signedKME(merged_datExpr, MEs)  # Rows = genes, Columns = MEs

# 3. Module number and color assignment
module_numbers <- merged_net$colors
module_colors <- labels2colors(module_numbers)

# 4. Map module numbers to actual ME column names
ME_names <- colnames(gene_kME)                     # e.g., "ME9", "ME8", etc.
ME_numbers <- as.numeric(sub("ME", "", ME_names))  # extract numeric part
ME_mapping <- setNames(ME_names, ME_numbers)      # module_number -> ME column

# 5. Vectorized assignment of kME for each gene
kME_values <- numeric(length(module_numbers))
for (i in seq_along(module_numbers)) {
  mod_num <- module_numbers[i]
  gene <- colnames(merged_datExpr)[i]
  
  if (mod_num %in% ME_numbers) {
    ME_col <- ME_mapping[as.character(mod_num)]
    kME_values[i] <- gene_kME[gene, ME_col]
  } else {
    # If the module number is missing (should not happen), take max absolute correlation
    kME_values[i] <- max(abs(gene_kME[gene, ]))
  }
}

# 6. Create output dataframe
kME_df <- data.frame(
  gene_name = colnames(merged_datExpr),
  module_number = module_numbers,
  module_color = module_colors,
  kME = kME_values,
  stringsAsFactors = FALSE
)

# 7. Sort by module number and descending kME
kME_df_sorted <- kME_df %>%
  arrange(module_number, desc(kME))

# 8. Save CSV
output_csv <- paste0("Merged_WGCNA_gene_kME.csv")
write.csv(kME_df_sorted, output_csv, row.names = FALSE)
cat(paste0(" -> Gene kME values saved to '", output_csv, "'\n"))


# --- Step 10: Save WGCNA Network Object ---
cat("\n### Saving WGCNA network object... ###\n")
# Save the essential components needed for downstream analysis (MEs)
# Note: Saving the full merged_net object can be large. MEs are often sufficient.
merged_MEs_output <- merged_net$MEs
save(merged_MEs_output, file = paste0(OUTPUT_PREFIX, "_MEs_object.RData"))
# Optionally save the full object if needed later:
# save(merged_net, file = paste0(OUTPUT_PREFIX, "_net_object.RData"))
cat(paste0(" -> WGCNA Module Eigengenes saved to '", OUTPUT_PREFIX, "_MEs_object.RData'.\n"))


cat("\n### WGCNA analysis for Merged dataset complete. ###\n")
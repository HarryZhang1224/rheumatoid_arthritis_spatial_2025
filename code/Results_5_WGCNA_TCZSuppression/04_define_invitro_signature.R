# ==================================================================
# R SCRIPT: All-vs-All Differential Expression Analysis (90 UP Signatures)
#
# Description:
# - Recreates the older 90-contrast IM049 workflow that matched the paper
#   scoring outputs.
# - Converts the one-hot IM049 metadata into a single Condition factor.
# - Fits one limma model across all 10 conditions.
# - Generates all 90 directional contrasts as "Condition_A - Condition_B".
# - Keeps the top 100 upregulated genes (logFC > 0, unadjusted P < 0.05)
#   for each contrast.
#
# Outputs:
# - 2026-03-19_audit_IM049_top100_up.csv
# ==================================================================

suppressPackageStartupMessages({
  library(limma)
})

EXPRESSION_FILE <- "IM049_data_proteincoding.csv"
METADATA_FILE <- "IM049_metadata.csv"
OUTPUT_FILE <- "2026-03-19_audit_IM049_top100_up.csv"

cat("### Loading IM049 expression and metadata... ###\n")

expr_data <- read.csv(EXPRESSION_FILE, row.names = 1, check.names = FALSE)
metadata_raw <- read.csv(METADATA_FILE, check.names = FALSE)
colnames(metadata_raw) <- trimws(colnames(metadata_raw))

condition_cols <- setdiff(colnames(metadata_raw), "SampleID")

# Convert one-hot metadata into a single condition label per sample.
metadata_conditions <- data.frame(
  SampleID_raw = metadata_raw$SampleID,
  Condition = NA_character_,
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(metadata_raw))) {
  hit_idx <- which(metadata_raw[i, condition_cols] == 1)

  if (length(hit_idx) != 1) {
    stop(
      paste0(
        "Each IM049 sample must map to exactly one condition. Problem at row ",
        i,
        " for sample ",
        metadata_raw$SampleID[i],
        "."
      )
    )
  }

  metadata_conditions$Condition[i] <- condition_cols[hit_idx]
}

# Align metadata to the expression matrix column names as read by R.
limma_metadata <- data.frame(
  SampleID = colnames(expr_data),
  stringsAsFactors = FALSE
)

limma_metadata$Condition <- metadata_conditions$Condition[
  match(limma_metadata$SampleID, make.names(metadata_conditions$SampleID_raw))
]

if (anyNA(limma_metadata$Condition)) {
  missing_samples <- limma_metadata$SampleID[is.na(limma_metadata$Condition)]
  stop(
    paste(
      "Failed to map IM049 expression columns back to metadata samples:",
      paste(missing_samples, collapse = ", ")
    )
  )
}

limma_metadata$Condition <- factor(limma_metadata$Condition)

cat("### Building limma design matrix... ###\n")

design <- model.matrix(~ 0 + Condition, data = limma_metadata)
colnames(design) <- levels(limma_metadata$Condition)

fit <- lmFit(expr_data, design)

condition_levels <- levels(limma_metadata$Condition)
condition_pairs <- combn(condition_levels, 2)

contrast_strings <- character()
for (i in seq_len(ncol(condition_pairs))) {
  pair <- condition_pairs[, i]
  contrast_strings <- c(
    contrast_strings,
    paste(pair[1], pair[2], sep = " - "),
    paste(pair[2], pair[1], sep = " - ")
  )
}

cat("### Generated ", length(contrast_strings), " directional contrasts. ###\n", sep = "")

contrast_matrix <- makeContrasts(contrasts = contrast_strings, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

all_up_genes <- list()

for (contrast_name in colnames(contrast_matrix)) {
  res_table <- topTable(fit2, coef = contrast_name, number = Inf, sort.by = "p")
  res_table_up <- res_table[res_table$logFC > 0 & res_table$P.Value < 0.05, , drop = FALSE]

  if (nrow(res_table_up) == 0) {
    next
  }

  res_table_top100 <- head(res_table_up, 100)
  res_table_top100$gene_name <- rownames(res_table_top100)
  res_table_top100$contrast <- contrast_name
  res_table_top100$direction <- "Upregulated"

  all_up_genes[[contrast_name]] <- res_table_top100
}

if (length(all_up_genes) == 0) {
  stop("No IM049 contrast produced any upregulated genes passing the filter.")
}

final_results <- do.call(rbind, all_up_genes)
final_results <- final_results[, c("contrast", "gene_name", "logFC", "P.Value", "adj.P.Val", "direction")]

write.csv(final_results, OUTPUT_FILE, row.names = FALSE)

cat("### IM049 all-vs-all signature generation complete. ###\n")
cat("Saved:", OUTPUT_FILE, "\n")
cat("Total rows:", nrow(final_results), "\n")

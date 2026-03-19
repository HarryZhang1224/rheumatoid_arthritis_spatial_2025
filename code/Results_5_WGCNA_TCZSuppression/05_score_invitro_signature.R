# ==================================================================
# R SCRIPT: Calculate Patient-Level Gene Set Scores and Deltas
#
# Description:
# - Loads the audit IM049 top-100 gene-set file produced by
#   04_define_invitro_signature.R.
# - Scores every IM049 gene set in paired TCZ patient bulk RNA-seq samples
#   as the average expression across genes in the set.
# - Calculates patient-level W16 - W0 deltas.
# - Prints the paper-relevant contrast so it can be checked quickly.
#
# Outputs:
# - 2026-03-19_audit_IM049_Top100_GeneSet_Deltas_TCZ.csv
# ==================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
})

IM049_DEG_FILE <- "2026-03-19_audit_IM049_top100_up.csv"
METADATA_FILE <- "R4RA_metadata.csv"
EXPR_FILE <- "TCZ_data_proteincoding.csv"
OUTPUT_FILE <- "2026-03-19_audit_IM049_Top100_GeneSet_Deltas_TCZ.csv"

PAPER_GENE_SET <- "M_CSFplusTGF_BplusFibs - M_CSFplusFibsplusTGF_BplusTCZ"

cat("### Loading audit IM049 gene sets... ###\n")

im049_data <- read.csv(IM049_DEG_FILE, check.names = FALSE)

im049_gene_sets <- im049_data %>%
  select(contrast, gene_name) %>%
  distinct() %>%
  split(.$contrast) %>%
  lapply(function(x) x$gene_name)

cat("Loaded ", length(im049_gene_sets), " IM049 gene sets.\n", sep = "")

cat("### Loading TCZ paired patient data... ###\n")

metadata_raw <- read.csv(METADATA_FILE, check.names = FALSE)
pheno_subset <- metadata_raw %>%
  filter(Drug == "TOC", Timepoint %in% c(0, 16))

paired_patient_ids <- pheno_subset %>%
  group_by(Patient_ID) %>%
  filter(n() == 2) %>%
  pull(Patient_ID) %>%
  unique()

pheno_subset_paired <- pheno_subset %>%
  filter(Patient_ID %in% paired_patient_ids)

datExpr <- read.csv(EXPR_FILE, check.names = FALSE) %>%
  column_to_rownames("gene_name") %>%
  as.matrix()

colnames(datExpr) <- sub("^X", "", colnames(datExpr))

datExpr_aligned <- datExpr[, colnames(datExpr) %in% pheno_subset_paired$Header, drop = FALSE]

if (ncol(datExpr_aligned) == 0) {
  stop("No paired TCZ patient samples were found in the expression matrix.")
}

score_list <- list()

for (set_name in names(im049_gene_sets)) {
  gene_set <- im049_gene_sets[[set_name]]
  genes_in_data <- intersect(gene_set, rownames(datExpr_aligned))

  if (length(genes_in_data) < 5) {
    next
  }

  subset_expr <- datExpr_aligned[genes_in_data, , drop = FALSE]
  avg_scores <- colMeans(subset_expr, na.rm = TRUE)

  score_list[[set_name]] <- data.frame(
    Header = names(avg_scores),
    Gene_Set = set_name,
    Avg_logCPM = avg_scores,
    stringsAsFactors = FALSE
  )
}

if (length(score_list) == 0) {
  stop("No IM049 gene sets had at least 5 overlapping genes in the TCZ data.")
}

tcz_scores <- bind_rows(score_list)

final_summary <- tcz_scores %>%
  left_join(
    pheno_subset_paired %>% select(Header, Patient_ID, Timepoint, Drug),
    by = "Header"
  ) %>%
  select(Patient_ID, Drug, Gene_Set, Timepoint, Avg_logCPM) %>%
  drop_na() %>%
  pivot_wider(
    names_from = Timepoint,
    values_from = Avg_logCPM,
    names_prefix = "W"
  ) %>%
  mutate(Delta_logCPM = W16 - W0) %>%
  drop_na(W0, W16) %>%
  arrange(Gene_Set, Patient_ID)

write.csv(final_summary, OUTPUT_FILE, row.names = FALSE)

cat("### TCZ patient scoring complete. ###\n")
cat("Saved:", OUTPUT_FILE, "\n")

paper_rows <- final_summary %>%
  filter(Gene_Set == PAPER_GENE_SET) %>%
  arrange(Patient_ID)

if (nrow(paper_rows) > 0) {
  cat("\n### Paper-relevant gene set preview ###\n")
  print(paper_rows[, c("Patient_ID", "W0", "W16", "Delta_logCPM")], row.names = FALSE)
}

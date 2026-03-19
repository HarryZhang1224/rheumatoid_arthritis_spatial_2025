# ==================================================================
# Analysis Code for "SPP1hi macrophages in fibrin niches..." Manuscript
# ==================================================================
#
# This repository contains the R scripts used for the Weighted Gene Co-expression
# Network Analysis (WGCNA), Module Eigengene (ME) analysis, functional annotation,
# and specialized gene signature scoring presented in the manuscript.
#
# The pipeline is organized into chronological steps, where the output of one script
# serves as the input for the next.
#
# ==================================================================

# --- OVERVIEW ---
# The analysis workflow involves three primary phases:
# 1. Data Preparation: Filtering raw expression data and constructing the universal WGCNA matrix.
# 2. Module Discovery & Annotation: Performing WGCNA and annotating modules via ORA.
# 3. Signature Validation: Scoring and statistically testing specific in vitro signatures
#    against paired in vivo patient changes (W16 vs W0).

## Requirements

# R Packages (required for running the scripts):
# NOTE: BiocManager is required to install WGCNA, limma, clusterProfiler, and ReactomePA.
# Example: BiocManager::install(c("limma", "WGCNA", "clusterProfiler", "ReactomePA"))

# R (version X.X.X)
# R Packages (with versions):
# * WGCNA (vX.X)
# * dplyr (vX.X)
# * tidyr (vX.X)
# * tibble (vX.X)
# * stringr (vX.X)
# * limma (vX.X)
# * ggplot2 (vX.X)
# * ggpubr (vX.X)
# * clusterProfiler (vX.X)
# * org.Hs.eg.db (vX.X)
# * ReactomePA (vX.X)
# * BiocManager (for installing Bioconductor packages)

## Input Data Files

# Ensure the following CSV files are in the same directory as the scripts,
# or modify the file paths within the scripts:
#
# * Merged_data_proteincoding.csv: Universal gene expression matrix (Genes x Samples)
# * TCZ_data_proteincoding.csv: Gene expression matrix for the TCZ dataset.
# * IM049_data_proteincoding.csv: Gene expression matrix for the in vitro experiment.
# * R4RA_metadata.csv: Metadata for the patient samples (TCZ and RTX).
# * IM049_metadata.csv: Metadata for the in vitro experiment samples.

# ==================================================================
## Analysis Pipeline
# ==================================================================

### Phase I: Data Preparation and Module Discovery

# **1. prepare_merged_protein_coding_data.R**
# ----------------------------------------
# **Purpose:** Data standardization. Filters all raw expression files to keep only
# **protein-coding genes** (using BiomaRt) and merges the individual TCZ and RTX
# datasets into one universal matrix.
# **Inputs:** *_data.csv (e.g., TCZ_data.csv, RTX_data.csv, IM049_data.csv)
# **Outputs:** Merged_data_proteincoding.csv

# **2. 01_run_wgcna_merged.R**
# ----------------------------
# **Purpose:** Weighted Gene Co-expression Network Analysis. Identifies robust co-expression
# patterns across the combined patient cohort.
# **Outputs:**
# * Merged_WGCNA_module_assignments.csv: Final gene-to-module assignments.
# * Merged_WGCNA_gene_kME.csv: Gene Module Membership (kME) values.

# **3. 02_analyze_module_changes.R**
# ---------------------------------
# **Purpose:** Module Eigengene (ME) analysis and clinical relevance (relevant to Fig 5A).
# * Calculates the **change in Module Eigengenes (ΔME = W16 - W0)** for the *Merged*
#   modules within each drug group (TCZ, RTX).
# * **Compares the ΔME distributions between TCZ and RTX groups** using an unpaired t-test.
# **Outputs:**
# * MergedModules_TCZ_RTX_Full_Delta_Comparison_Final.csv
# * MergedModules_in_TCZ_paired_analysis.csv
# * MergedModules_in_RTX_paired_analysis.csv

### Phase II: Functional Annotation

# **4. 03_run_ora_merged_modules.R**
# ----------------------------------
# **Purpose:** Functional annotation of WGCNA modules.
# * Performs Over-Representation Analysis (ORA) on the **hub genes** (kME > 0.5 or kME < -0.5)
#   for each *Merged* WGCNA module.
# * **Saves all significant results (Adj. P < 0.05)** for each database.
# **Inputs:** Merged_WGCNA_gene_kME.csv
# **Outputs:** ORA_Significant_GOBP_ModuleHubs.csv, ORA_Significant_KEGG_ModuleHubs.csv, ORA_Significant_Reactome_ModuleHubs.csv

### Phase III: Signature Validation and Paired Testing

# **5. 04_define_invitro_signature.R**
# ------------------------------------
# **Purpose:** Define the IM049 in vitro top-100 signatures used for patient scoring.
# * Performs limma differential expression across all 10 IM049 in vitro conditions.
# * Converts the one-hot metadata into a single condition factor and fits one model.
# * Generates all 90 directional contrasts and keeps the top 100 upregulated genes
#   (logFC > 0, unadjusted P < 0.05) for each contrast.
# **Outputs:** 2026-03-19_audit_IM049_top100_up.csv

# **6. 05_score_invitro_signature.R**
# ------------------------------------
# **Purpose:** Score all IM049 top-100 signatures in paired TCZ patient samples.
# * Loads the IM049 contrast gene sets from script 04.
# * Calculates the **average logCPM score** for each gene set in each paired TOC sample.
# * Reshapes the results to paired Week 0 / Week 16 format and calculates W16 - W0 deltas.
# * Prints the paper-relevant contrast (`M_CSFplusTGF_BplusFibs - M_CSFplusFibsplusTGF_BplusTCZ`)
#   as a quick audit check.
# **Outputs:** 2026-03-19_audit_IM049_Top100_GeneSet_Deltas_TCZ.csv

# **7. wilcoxon_paired_analysis_90_sets.R**
# ------------------------------------------
# **Purpose:** Comprehensive statistical test for all 90 signatures. (Note: This script is typically run for exploratory analysis to generate the full list of significant pathways.)
# * Loads all 90 IM049 UP signatures.
# * Runs the final **Wilcoxon Matched-Pairs Signed-Rank Test** on the Deltas of all 90
#   gene sets and applies BH correction.
# **Outputs:** Paired_DEGs_Summary_TCZ_Wilcoxon_BH_90_Sets.csv

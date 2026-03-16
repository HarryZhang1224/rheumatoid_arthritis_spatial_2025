# ==================================================================
# SCRIPT: ORA for Merged WGCNA Module Hub Genes
#
# Description:
# - Performs Over-Representation Analysis (ORA) for GO:BP, KEGG, and
#   Reactome pathways.
# - Analyzes hub genes (kME > 0.5 or kME < -0.5) for each module
#   identified in the *Merged* WGCNA analysis.
# - Saves all significant results (Adjusted P < 0.05) per database.
#
# Inputs:
# - Merged_WGCNA_gene_kME.csv: Gene module assignments and kME values.
#
# Outputs:
# - ORA_Significant_GOBP_ModuleHubs.csv: Significant GO:BP enrichments.
# - ORA_Significant_KEGG_ModuleHubs.csv: Significant KEGG enrichments.
# - ORA_Significant_Reactome_ModuleHubs.csv: Significant Reactome enrichments.
# ==================================================================


# ==================================================================
# SCRIPT: ORA for Merged WGCNA Module Hub Genes
# ==================================================================

# --- Step 1: Load Libraries ---
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)

# --- Step 2: Define File Paths and Parameters ---
KME_FILE <- "Merged_WGCNA_gene_kME.csv"
DATASET_NAME <- "Merged"
OUTPUT_PREFIX <- "ORA_Significant"

# --- Step 3: Load kME Data and Define Gene Universe ---
cat("### Loading WGCNA kME data and defining gene universe... ###\n")
kme_data <- read.csv(KME_FILE)
universe_genes_symbol <- unique(kme_data$gene_name)
cat(paste0(" -> Universe size (Symbols): ", length(universe_genes_symbol), " genes.\n"))

# --- Step 4: Convert Gene Symbols to Entrez IDs ---
cat(" -> Converting gene symbols to Entrez IDs...\n")
universe_ids <- tryCatch({
  bitr(universe_genes_symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
}, warning = function(w){
  message("Warning during universe ID conversion: ", conditionMessage(w))
  suppressWarnings(bitr(universe_genes_symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db))
}, error = function(e){
  stop("Fatal Error during universe ID conversion: ", conditionMessage(e))
})
universe_entrez <- unique(universe_ids$ENTREZID[!is.na(universe_ids$ENTREZID)])
cat(paste0(" -> Universe size (Valid Unique Entrez IDs): ", length(universe_entrez), " genes.\n"))

# --- Step 5: Prepare Module Hub Gene Lists (Entrez IDs) ---
cat("### Preparing module hub gene lists (Entrez IDs)... ###\n")
module_hub_entrez_lists <- list()
modules_to_analyze <- unique(kme_data$module_color[kme_data$module_color != "grey"])

for (mod_color in modules_to_analyze) {
  module_subset <- kme_data %>% dplyr::filter(module_color == mod_color)
  
  # Positive hubs
  genes_symbol_pos <- module_subset %>% dplyr::filter(kME > 0.5) %>% dplyr::pull(gene_name)
  if (length(genes_symbol_pos) >= 5) {
    ids_pos <- tryCatch({
      bitr(genes_symbol_pos, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID
    }, warning = function(w) {
      message("Warning converting Pos hubs for module ", mod_color, ": ", conditionMessage(w))
      suppressWarnings(bitr(genes_symbol_pos, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID)
    }, error = function(e) { message("Error converting Pos hubs for module ", mod_color); NULL })
    
    ids_pos_clean <- unique(ids_pos[!is.na(ids_pos)])
    if (length(ids_pos_clean) >= 5) {
      hub_name_pos <- paste0("ME", mod_color, "_Pos")
      module_hub_entrez_lists[[hub_name_pos]] <- ids_pos_clean
      cat(paste0("  - ", hub_name_pos, ": ", length(ids_pos_clean), " Entrez IDs.\n"))
    }
  }
  
  # Negative hubs
  genes_symbol_neg <- module_subset %>% dplyr::filter(kME < -0.5) %>% dplyr::pull(gene_name)
  if (length(genes_symbol_neg) >= 5) {
    ids_neg <- tryCatch({
      bitr(genes_symbol_neg, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID
    }, warning = function(w) {
      message("Warning converting Neg hubs for module ", mod_color, ": ", conditionMessage(w))
      suppressWarnings(bitr(genes_symbol_neg, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID)
    }, error = function(e) { message("Error converting Neg hubs for module ", mod_color); NULL })
    
    ids_neg_clean <- unique(ids_neg[!is.na(ids_neg)])
    if (length(ids_neg_clean) >= 5) {
      hub_name_neg <- paste0("ME", mod_color, "_Neg")
      module_hub_entrez_lists[[hub_name_neg]] <- ids_neg_clean
      cat(paste0("  - ", hub_name_neg, ": ", length(ids_neg_clean), " Entrez IDs.\n"))
    }
  }
}

if (length(module_hub_entrez_lists) == 0) stop("ERROR: No module hub lists with >=5 Entrez IDs.")

# --- Step 6: Define ORA Function ---
run_ora <- function(gene_list, hub_name, universe_list, analysis_type) {
  cat(paste0("    * Running ", analysis_type, " for ", hub_name, "... "))
  enrich_result <- NULL
  
  if (length(gene_list) < 5 || length(universe_list) < 10) return(NULL)
  gene_list_in_universe <- intersect(gene_list, universe_list)
  if (length(gene_list_in_universe) < 5) return(NULL)
  
  tryCatch({
    if (analysis_type == "GO:BP") {
      enrich_result <- enrichGO(
        gene = gene_list_in_universe,
        universe = universe_list,
        OrgDb = org.Hs.eg.db,
        keyType = "ENTREZID",
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = 1,
        qvalueCutoff = 1
      )
    } else if (analysis_type == "KEGG") {
      enrich_result <- enrichKEGG(
        gene = gene_list_in_universe,
        universe = universe_list,
        organism = 'hsa',
        pAdjustMethod = "BH",
        pvalueCutoff = 1,
        qvalueCutoff = 1,
        minGSSize = 5
      )
    } else if (analysis_type == "Reactome") {
      enrich_result <- ReactomePA::enrichPathway(
        gene = gene_list_in_universe,
        universe = universe_list,
        organism = "human",
        pAdjustMethod = "BH",
        pvalueCutoff = 1,
        qvalueCutoff = 1,
        minGSSize = 5
      )
    }
  }, error = function(e) { cat("ERROR:", conditionMessage(e), "\n") })
  
  if (!is.null(enrich_result) && nrow(enrich_result@result) > 0) {
    cat("Done. ")
    if (analysis_type != "KEGG") {
      tryCatch({
        enrich_result <- setReadable(enrich_result, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
        cat("(Readable) ")
      }, warning = function(w){}, error = function(e){})
    }
    cat("\n")
    
    results_df <- as.data.frame(enrich_result) %>%
      dplyr::mutate(
        Module_Hub_Set = hub_name,
        Database = analysis_type,
        Adjusted_P_value = p.adjust,
        Overlap_Genes = if("geneID" %in% names(.)) gsub("/", ", ", geneID) else NA_character_
      ) %>%
      dplyr::select(
        Module_Hub_Set,
        Database,
        ID,
        Description,
        Adjusted_P_value,
        Overlap_Genes,
        Count
      ) %>%
      dplyr::filter(!is.na(Adjusted_P_value))
    
    return(results_df)
  } else {
    return(NULL)
  }
}

# --- Step 7: Run ORA and Combine Results ---
cat("\n### Running ORA across databases and module lists... ###\n")
databases <- c("GO:BP", "KEGG", "Reactome")
all_database_results <- list()

for (db in databases) {
  cat(paste0(" -> Database: ", db, "\n"))
  db_results_list <- list()
  
  for (hub_name in names(module_hub_entrez_lists)) {
    gene_list_entrez <- module_hub_entrez_lists[[hub_name]]
    if (length(gene_list_entrez) < 5) next
    ora_df <- run_ora(gene_list_entrez, hub_name, universe_entrez, db)
    if (!is.null(ora_df) && nrow(ora_df) > 0) db_results_list[[hub_name]] <- ora_df
  }
  
  if (length(db_results_list) > 0) all_database_results[[db]] <- dplyr::bind_rows(db_results_list)
}

# --- Step 8: Save Significant Results ---
cat("\n### Filtering significant results (Adj. P < 0.05) and saving CSVs... ###\n")
for (db in names(all_database_results)) {
  db_combined_df <- all_database_results[[db]] %>%
    dplyr::filter(Adjusted_P_value < 0.05) %>%
    dplyr::arrange(Module_Hub_Set, Adjusted_P_value)
  
  if (nrow(db_combined_df) > 0) {
    db_filename_part <- gsub("[: ]", "", db)
    output_filename <- paste0(OUTPUT_PREFIX, "_", db_filename_part, "_ModuleHubs.csv")
    write.csv(db_combined_df, output_filename, row.names = FALSE)
    cat(paste0(" -> Saved ", nrow(db_combined_df), " significant results to '", output_filename, "'\n"))
  } else {
    cat(paste0(" -> No significant results for database: ", db, "\n"))
  }
}

cat("\n### All ORA analyses complete. ###\n")


# -------------------------------
# POST-PROCESS ORA CSVs TO ADD MODULE NUMBERS
# -------------------------------

# Load necessary library
library(stringr) # for str_extract

# List of ORA CSVs to update
ora_files <- c(
  "ORA_Significant_GOBP_ModuleHubs.csv",
  "ORA_Significant_KEGG_ModuleHubs.csv",
  "ORA_Significant_Reactome_ModuleHubs.csv"
)

# Load kME data for module color -> number mapping
kme_data <- read.csv("Merged_WGCNA_gene_kME.csv")
color_to_number <- setNames(
  kme_data$module_number,
  kme_data$module_color
)
color_to_number <- color_to_number[!duplicated(names(color_to_number))]

# Function to process a single ORA CSV
add_module_number <- function(ora_file, color_map) {
  if (!file.exists(ora_file)) {
    message("File not found: ", ora_file)
    return(NULL)
  }
  
  ora_df <- read.csv(ora_file)
  
  # Extract module color from Module_Hub_Set
  ora_df$Module_Color <- str_extract(ora_df$Module_Hub_Set, "(?<=ME)[A-Za-z]+")
  
  # Map color to module number
  ora_df$Module_Number <- color_map[ora_df$Module_Color]
  
  # Reorder columns to put number after color
  ora_df <- ora_df[, c(
    "Module_Hub_Set",
    "Module_Color",
    "Module_Number",
    setdiff(names(ora_df), c("Module_Hub_Set", "Module_Color", "Module_Number"))
  )]
  
  # Save back to CSV (overwrite original)
  write.csv(ora_df, ora_file, row.names = FALSE)
  message("Updated ORA CSV with module numbers saved to '", ora_file, "'")
}

# Loop through all ORA CSVs
lapply(ora_files, add_module_number, color_map = color_to_number)

# ==================================================================
# SCRIPT: Analyze Paired Changes in Merged WGCNA Module Eigengenes
#
# Description:
# - Calculates the change in Module Eigengenes (ΔME = W16 - W0) for
#   modules identified in the *Merged* WGCNA analysis.
# - Performs calculations separately within TCZ and RTX patient groups.
# - Compares the distribution of ΔME scores between TCZ and RTX using
#   an unpaired t-test (relevant for Fig 5A).
# - Saves summary statistics and generates magnitude plots.
#
# Inputs:
# - Merged_WGCNA_MEs_object.RData: Contains the 'merged_MEs_output' object from script 01.
# - R4RA_metadata.csv: Metadata for the patient samples.
#
# Outputs:
# - MergedModules_in_TCZ_paired_analysis.csv: Paired ΔME stats within TCZ group.
# - MergedModules_in_RTX_paired_analysis.csv: Paired ΔME stats within RTX group.
# - MergedModules_TCZ_RTX_Full_Delta_Comparison_Final.csv: Comparison of ΔME between TCZ and RTX.
# - MergedModules_in_TCZ_magnitude_plot.png: Plot of ΔME magnitude for TCZ group.
# - MergedModules_in_RTX_magnitude_plot.png: Plot of ΔME magnitude for RTX group.
# ==================================================================

# --- Step 1: Load Libraries ---
# if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
# if (!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr")
# if (!requireNamespace("tibble", quietly = TRUE)) install.packages("tibble")
# if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
# if (!requireNamespace("gtools", quietly = TRUE)) install.packages("gtools") # For mixedsort

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(gtools) # For mixed sorting of module colors/numbers

# --- Step 2: Define File Paths and Load Data ---
WGCNA_ME_FILE <- "Merged_WGCNA_MEs_object.RData" # File containing Module Eigengenes
METADATA_FILE <- "R4RA_metadata.csv"
OUTPUT_PREFIX <- "MergedModules"

cat("### Loading WGCNA Module Eigengenes and metadata... ###\n")

# Load the WGCNA Module Eigengenes object ('merged_MEs_output')
load(WGCNA_ME_FILE) # This should load the 'merged_MEs_output' object

# Convert MEs to a dataframe with Header column if necessary
merged_MEs_df <- merged_MEs_output %>%
  as.data.frame() %>%
  rownames_to_column("Header")

# Load metadata
R4RA_metadata <- read.csv(METADATA_FILE)

# --- Step 3: Define Reusable Analysis Function ---
perform_paired_analysis_by_drug <- function(drug_code, drug_name, all_MEs_df, metadata) {
  
  cat(paste0("\n### Analyzing Merged Modules in ", drug_name, " Patients ###\n"))
  
  # --- Filter Data for the Specific Drug Group ---
  ME_data_subset <- metadata %>%
    filter(Drug == drug_code) %>%
    select(Header, Patient_ID, Timepoint) %>%
    inner_join(all_MEs_df, by = "Header")
  
  # Calculate the change (ΔME) for each patient in this group
  # Reshape to get W0 and W16 side-by-side
  delta_MEs_subset <- ME_data_subset %>%
    select(-Header) %>% # Drop original sample ID
    pivot_longer(cols = starts_with("ME"), names_to = "Module", values_to = "ME_value") %>%
    pivot_wider(names_from = Timepoint, values_from = ME_value, names_prefix = "W") %>%
    # Keep only patients with both time points
    filter(!is.na(W0) & !is.na(W16)) %>%
    mutate(delta_ME = W16 - W0)
  
  if(nrow(delta_MEs_subset) == 0) {
    cat(paste0(" -> No paired samples found for ", drug_name, ". Skipping paired analysis.\n"))
    return(NULL) # Return NULL if no paired data
  }
  
  
  # --- Perform Paired T-Test on Each Module (W16 vs W0 within this drug group) ---
  cat(" -> Performing paired t-tests (W16 vs W0)...\n")
  module_stats_subset <- delta_MEs_subset %>%
    group_by(Module) %>%
    summarize(
      N_Paired = sum(!is.na(delta_ME)), # Count valid pairs
      if (N_Paired >= 2) {
        # Perform one-sample t-test on the difference (equivalent to paired t-test)
        test_result <- t.test(delta_ME, mu = 0)
        tibble(
          mean_delta_ME = test_result$estimate,
          t_statistic = test_result$statistic,
          p_value = test_result$p.value
        )
      } else {
        # Return NA if fewer than 2 pairs
        tibble(
          mean_delta_ME = mean(delta_ME, na.rm = TRUE), # Will be NaN if N < 1
          t_statistic = NA_real_,
          p_value = NA_real_
        )
      },
      .groups = 'drop'
    ) %>%
    filter(!is.na(p_value)) %>% # Remove modules where test couldn't run
    mutate(padj = p.adjust(p_value, method = "BH")) %>% # Adjust p-values across modules
    mutate(
      module_color = gsub("ME", "", Module),
      # Extract number, handle grey (0) correctly for sorting
      module_number = as.numeric(gsub("[^0-9]", "", module_color)), # Extract numeric part
      module_number = ifelse(module_color == "grey", 0, module_number) # Assign 0 to grey explicitly
    ) %>%
    arrange(module_number) # Sort by module number
  
  # --- Save CSV Output ---
  csv_filename <- paste0(OUTPUT_PREFIX, "_in_", drug_name, "_paired_analysis.csv")
  write.csv(module_stats_subset, csv_filename, row.names = FALSE)
  cat(paste(" -> Paired analysis results saved to '", csv_filename, "'\n"))
  
  # --- Generate and Save Magnitude Plot ---
  plot_data <- module_stats_subset %>%
    mutate(
      Direction = ifelse(mean_delta_ME > 0, "Upregulated", "Downregulated"),
      Magnitude = abs(mean_delta_ME),
      significance = case_when(padj < 0.05 ~ "*", p_value < 0.05 ~ "+", TRUE ~ ""),
      Label_Vjust = ifelse(mean_delta_ME > 0, -0.5, 1.5) # Adjust label position based on bar direction (relative to magnitude)
    ) %>%
    # Exclude the grey module (ME0) from plotting if desired
    filter(module_color != "grey")
  
  if (nrow(plot_data) > 0) {
    # Define custom colors for plot
    CUSTOM_COLORS <- c("Upregulated" = "#B80000", "Downregulated" = "#0072B2") # Red/Blue
    
    plot_object <- ggplot(plot_data, aes(x = reorder(module_color, module_number), y = Magnitude, fill = Direction)) +
      geom_bar(stat = "identity", color = "black", width = 0.8) +
      scale_fill_manual(values = CUSTOM_COLORS) +
      # Use Magnitude for y position, adjust vjust
      geom_text(aes(label = significance, y = Magnitude), vjust = plot_data$Label_Vjust, size = 6) +
      geom_hline(yintercept = 0) +
      labs(
        title = paste("Magnitude of Merged Module Change (", drug_name, " Patients: W16 vs W0)"),
        subtitle = "Color: Direction (Red=UP, Blue=DOWN). *: Adj P<0.05, +: P<0.05",
        x = "WGCNA Module (from Merged Analysis)",
        y = "Magnitude of Mean Eigengene Change (|ΔME|)"
      ) +
      theme_classic(base_size = 14) +
      theme(
        axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.title = element_blank()
      )
    
    plot_filename <- paste0(OUTPUT_PREFIX, "_in_", drug_name, "_magnitude_plot.png")
    ggsave(plot_filename, plot = plot_object, width = 14, height = 7, dpi = 300)
    cat(paste(" -> Magnitude plot saved to '", plot_filename, "'\n"))
  } else {
    cat(" -> Skipping magnitude plot: No non-grey modules with results.\n")
  }
  
  
  # Return the delta MEs for the comparison step
  return(delta_MEs_subset %>% select(Patient_ID, Module, delta_ME) %>% mutate(Drug = drug_code))
}

# --- Step 4: Run the Paired Analysis for Both Drug Groups ---
tcz_deltas_long <- perform_paired_analysis_by_drug(drug_code = "TOC", drug_name = "TCZ", all_MEs_df = merged_MEs_df, metadata = R4RA_metadata)
rtx_deltas_long <- perform_paired_analysis_by_drug(drug_code = "RTX", drug_name = "RTX", all_MEs_df = merged_MEs_df, metadata = R4RA_metadata)

# --- Step 5: Compare Delta MEs between TCZ and RTX (Unpaired T-test) ---
# Check if both analyses returned results
if (!is.null(tcz_deltas_long) && !is.null(rtx_deltas_long)) {
  
  cat("\n### Comparing Delta MEs between TCZ and RTX groups... ###\n")
  
  # Combine individual patient deltas from both groups
  combined_deltas_long <- bind_rows(tcz_deltas_long, rtx_deltas_long) %>%
    filter(!is.na(delta_ME)) %>% # Ensure only valid deltas are used
    mutate(Drug = factor(Drug)) # Convert Drug to factor for t.test
  
  # Perform unpaired t-test for each module, comparing delta_ME between Drug groups
  unpaired_t_test_results <- combined_deltas_long %>%
    group_by(Module) %>%
    summarise(
      N_TCZ = sum(Drug == "TOC"),
      N_RTX = sum(Drug == "RTX"),
      # Run t.test only if both groups have at least 2 samples
      if (N_TCZ >= 2 & N_RTX >= 2) {
        # Welch's t-test (var.equal = FALSE) is generally safer
        # Use tryCatch in case t.test fails for unforeseen reasons (e.g., zero variance within a group)
        test_result <- tryCatch({
          t.test(delta_ME ~ Drug, data = cur_data(), var.equal = FALSE)
        }, error = function(e) NULL)
        
        if(!is.null(test_result)){
          tibble(
            T_Statistic_Unpaired = test_result$statistic,
            P_value_Unpaired = test_result$p.value
          )
        } else {
          tibble( T_Statistic_Unpaired = NA_real_, P_value_Unpaired = NA_real_)
        }
        
      } else {
        # Return NA if insufficient data
        tibble(
          T_Statistic_Unpaired = NA_real_,
          P_value_Unpaired = NA_real_
        )
      },
      .groups = 'drop'
    ) %>%
    filter(!is.na(P_value_Unpaired)) # Remove modules where test couldn't run
  
  # --- Step 6: Consolidate All Results and Save ---
  if (nrow(unpaired_t_test_results) > 0) {
    cat(" -> Consolidating all comparison results...\n")
    
    # Load the paired results again for merging
    tcz_paired_summary <- read.csv(paste0(OUTPUT_PREFIX, "_in_TCZ_paired_analysis.csv")) %>%
      select(Module, Mean_Delta_TCZ = mean_delta_ME, P_value_TCZ = p_value, Adj_P_value_TCZ = padj)
    
    rtx_paired_summary <- read.csv(paste0(OUTPUT_PREFIX, "_in_RTX_paired_analysis.csv")) %>%
      select(Module, Mean_Delta_RTX = mean_delta_ME, P_value_RTX = p_value, Adj_P_value_RTX = padj)
    
    # Join paired summaries with the unpaired comparison results
    final_comparison_summary_df <- unpaired_t_test_results %>%
      inner_join(tcz_paired_summary, by = "Module") %>%
      inner_join(rtx_paired_summary, by = "Module") %>%
      # Calculate the difference between the mean deltas
      mutate(Difference_TCZ_RTX = Mean_Delta_TCZ - Mean_Delta_RTX) %>%
      # Apply BH correction to the unpaired P-values
      mutate(Adj_P_value_Unpaired = p.adjust(P_value_Unpaired, method = "BH")) %>%
      # Add module color/number back for sorting/readability
      mutate(
        module_color = gsub("ME", "", Module),
        module_number = as.numeric(gsub("[^0-9]", "", module_color)), # Extract numeric part
        module_number = ifelse(module_color == "grey", 0, module_number) # Assign 0 to grey
      ) %>%
      # Order columns and sort rows
      select(
        Module, module_color, module_number,
        N_TCZ, N_RTX,
        Mean_Delta_TCZ, P_value_TCZ, Adj_P_value_TCZ,
        Mean_Delta_RTX, P_value_RTX, Adj_P_value_RTX,
        Difference_TCZ_RTX,
        T_Statistic_Unpaired, P_value_Unpaired, Adj_P_value_Unpaired
      ) %>%
      arrange(Adj_P_value_Unpaired, module_number) # Sort primarily by significance of difference
    
    # Save the final consolidated table
    comparison_csv_file <- paste0(OUTPUT_PREFIX, "_TCZ_RTX_Full_Delta_Comparison_Final.csv")
    write.csv(final_comparison_summary_df, comparison_csv_file, row.names = FALSE)
    
    cat(paste0(" -> Comprehensive comparison saved to '", comparison_csv_file, "'\n"))
    cat("\n--- Top Modules Differentially Changed Between TCZ and RTX ---\n")
    print(head(final_comparison_summary_df, 10))
  } else {
    cat(" -> No modules showed significant differences between TCZ and RTX delta MEs.\n")
  }
  
} else {
  cat("\n### Skipping comparison between TCZ and RTX delta MEs due to missing results from one or both groups. ###\n")
}


cat("\n### Module change analysis complete. ###\n")

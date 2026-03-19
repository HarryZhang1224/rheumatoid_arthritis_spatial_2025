# Results 5 Analysis Code

This directory contains the R scripts used for the Results 5 analyses in the manuscript.

## Script order

1. `00_prepare_data.R`
2. `01_run_wgcna_merged.R`
3. `02_analyze_module_changes.R`
4. `03_run_ora_merged_modules.R`
5. `04_define_invitro_signature.R`
6. `05_score_invitro_signature.R`

## What each script does

### `00_prepare_data.R`
Prepares the protein-coding expression matrices used in downstream analyses.

### `01_run_wgcna_merged.R`
Runs WGCNA on the merged TCZ and RTX synovial bulk RNA-seq dataset and exports module assignments, eigengenes, and kME values.

### `02_analyze_module_changes.R`
Calculates paired module eigengene changes within treatment groups and compares TCZ and RTX module deltas.

### `03_run_ora_merged_modules.R`
Performs over-representation analysis on WGCNA module hub genes.

### `04_define_invitro_signature.R`
Builds the IM049 in vitro gene-set collection used for Fig. 5E-style scoring. The script fits a limma model across all 10 IM049 conditions, generates all 90 directional contrasts, and saves the top 100 upregulated genes for each contrast.

Output:
- `2026-03-19_audit_IM049_top100_up.csv`

### `05_score_invitro_signature.R`
Scores the IM049 contrast gene sets in paired TCZ patient bulk RNA-seq samples by average logCPM and calculates `Week16 - Week0` deltas. The paper-relevant contrast is:

- `M_CSFplusTGF_BplusFibs - M_CSFplusFibsplusTGF_BplusTCZ`

Output:
- `2026-03-19_audit_IM049_Top100_GeneSet_Deltas_TCZ.csv`

## Expected inputs

These scripts were written to run from the parent analysis directory, where the input and output files live. In the manuscript working tree, that parent directory is:

- `MZ_Results_5_Analysis/`

Key expected input files include:

- `IM049_data_proteincoding.csv`
- `IM049_metadata.csv`
- `Merged_data_proteincoding.csv`
- `R4RA_metadata.csv`
- `TCZ_data_proteincoding.csv`
- `RTX_data_proteincoding.csv`

## Notes

- The active Fig. 5E pipeline in this directory is the paper-matching version.
- Older in vitro signature scripts that did not reproduce the final paper values were archived outside this code directory.

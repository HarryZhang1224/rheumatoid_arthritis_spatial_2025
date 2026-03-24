IF signal quantification script for submission support

Primary script
- quantify_if_signal.py

What it does
- Quantifies pSTAT3 signal in DAPI-anchored events from a pSTAT3 panel.
- Quantifies the proportion of CD14+OPN+ events assigned to fibrin-only, collagen-only, both, or neither ECM regions from a fibrin/collagen panel.

Required inputs
- One pSTAT3-panel figure-parameter CSV.
- One fibrin/collagen-panel figure-parameter CSV.
- A directory containing the matching OME-TIFF files.

Typical command
python quantify_if_signal.py \
  --pstat3-csv /path/to/pSTAT3_panel.figure_parameters.csv \
  --ecm-csv /path/to/fibrin_collagen_panel.figure_parameters.csv \
  --ome-search-root /path/to/raw_ome_tiffs \
  --out-dir /path/to/output_dir

Optional mask overrides
- --pstat3-mask-threshold 0.14
- --ecm-mask-threshold 0.14

Main outputs
- if_signal_quantification_summary.csv
- pstat3_seed_calls.csv
- ecm_seed_calls.csv
- qc/*.mask_qc.png

Notes for manuscript support
- The script is generic and does not contain patient-specific identifiers in its logic.
- It is intended to document and reproduce the two IF quantifications reported in the manuscript.
- See README.md for a GitHub-facing overview and installation instructions.

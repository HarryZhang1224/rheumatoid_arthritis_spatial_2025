# IF Signal Quantification Submission Package

This repository-ready package contains the manuscript-facing code used for the two immunofluorescence quantifications reported in the paper:

1. pSTAT3 binning in the pSTAT3 panel
2. ECM compartment assignment of CD14+OPN+ events in the fibrin/collagen panel

## Contents

- `quantify_if_signal.py`
  Main analysis script.
- `requirements.txt`
  Tested Python package versions for this package.
- `PACKAGE_OVERVIEW.txt`
  Short package overview.
- `QUANTIFY_IF_SIGNAL_README.txt`
  Inputs, outputs, and method summary.
- `RUN_QUANTIFY_IF_SIGNAL_EXAMPLE.txt`
  Example command-line invocation.
- `regression_validation/`
  Regression evidence showing that the fixed script preserves reference-case outputs while failing loudly on the reviewed edge cases.

## Installation

Create and activate a Python environment, then install dependencies:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Typical usage

```bash
python quantify_if_signal.py \
  --pstat3-csv /path/to/pSTAT3_panel.figure_parameters.csv \
  --ecm-csv /path/to/fibrin_collagen_panel.figure_parameters.csv \
  --ome-search-root /path/to/raw_ome_tiffs \
  --out-dir /path/to/output_dir
```

Optional mask overrides:

```bash
python quantify_if_signal.py \
  --pstat3-csv /path/to/pSTAT3_panel.figure_parameters.csv \
  --ecm-csv /path/to/fibrin_collagen_panel.figure_parameters.csv \
  --ome-search-root /path/to/raw_ome_tiffs \
  --out-dir /path/to/output_dir \
  --pstat3-mask-threshold 0.14 \
  --ecm-mask-threshold 0.14
```

## Main outputs

- `if_signal_quantification_summary.csv`
- `pstat3_seed_calls.csv`
- `ecm_seed_calls.csv`
- `qc/*.mask_qc.png`

## Regression validation

See `regression_validation/REGRESSION_SUMMARY.txt`.

The packaged regression results show:

- identical outputs for the fixed script versus the pre-fix script on the reference case
- loud failure for bad marker CSVs
- loud failure for duplicate OME filename matches
- successful execution on single-level OME-TIFF inputs

## Notes

- The script logic is generic and does not encode patient-specific identifiers.
- The code package is intended to document and reproduce the manuscript IF quantifications rather than serve as a general-purpose image-analysis framework.

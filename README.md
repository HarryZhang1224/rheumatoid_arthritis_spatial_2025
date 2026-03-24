# SPP1hi macrophages in fibrin niches promote hyperplastic tissue remodeling in rheumatoid arthritis synovium


---
## Overview
This repository contains code to reproduce the main data processing, analysis, and figures for the manuscript *SPP1hi macrophages in fibrin niches promote hyperplastic tissue remodeling in rheumatoid arthritis synovium*.

1. **Co-culture scRNA-seq preprocessing** in [`code/scRNA/preprocessing_0`](./code/scRNA/preprocessing_0), including quality control, doublet detection, fibroblast/macrophage separation, and Symphony-based reference mapping to the Zhang et al. 2023 AMP atlas.
2. **Xenium in situ analysis** in [`code/xenium`](./code/xenium), including Baysor-based segmentation, Seurat-based quality control and integration, cell type niche identification, Niche-DE analysis, and cell type colocalization quotient (CLQ) analysis.
3. **Immunofluorescence analysis** in [`code/immunofluorescence_analysis`](./code/immunofluorescence_analysis), including the Python-based pipeline for pSTAT3 binning and ECM compartment assignment of CD14+OPN+ events from multiplex OME-TIFF images, together with regression-validation materials.
4. **Results 5 WGCNA and TCZ suppression analysis** in [`code/Results_5_WGCNA_TCZSuppression`](./code/Results_5_WGCNA_TCZSuppression), including merged WGCNA, paired module-change analysis, over-representation analysis, IM049 in vitro signature definition, and TCZ patient signature scoring.
5. Jupyter notebooks to reproduce **main manuscript figures** are in [`code/figures`](./code/figures).
---

## Data Access

Processed and raw data are available on [Synapse](https://www.synapse.org/Synapse:syn73710965/datasets/).

---

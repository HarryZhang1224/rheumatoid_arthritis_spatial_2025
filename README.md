# SPP1hi macrophages in fibrin niches promote hyperplastic tissue remodeling in rheumatoid arthritis synovium


---
## Overview
This repository contains code to reproduce the main data processing, analysis, and figures for the manuscript *SPP1hi macrophages in fibrin niches promote hyperplastic tissue remodeling in rheumatoid arthritis synovium*.

1. **Co-culture scRNA-seq preprocessing** in [`code/scRNA/preprocessing_0`](./code/scRNA/preprocessing_0), including quality control, doublet detection, fibroblast/macrophage separation, and Symphony-based reference mapping to the Zhang et al. 2023 AMP atlas.
2. **Xenium in situ analysis** in [`code/xenium`](./code/xenium), including Baysor-based segmentation, Seurat-based quality control and integration, cell type niche identification, Niche-DE analysis, and cell type colocalization quotient (CLQ) analysis.
3. **WGCNA** analysis in [`code/WGCNA_rheumatoid_arthritis`](./code/WGCNA_rheumatoid_arthritis).
3. Jupyter notebooks to reproduce **main manucript figures** are in [`code/figures`](./code/figures).
---

## Data Access

Processed and raw data are available on [Synapse](https://www.synapse.org/Synapse:syn73710965/datasets/).

---

# Xenium In Situ data quality control, integration, and subclustering(with Seurat v5 )

---
## Installation
Refer to [Seurat GitHub](https://github.com/satijalab/seurat).

---
## Pipeline Steps

1. **QC** (`QC_0`)
- The scripts here runs quality control on the Seurat objects created from Baysor outputs.

2. **Integration and subclustering** (`integration_clustering_1`)
- The scripts here runs atomic sketch integration, SingleR cell type annotations, and subclustering.
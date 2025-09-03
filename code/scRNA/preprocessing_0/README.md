# Preprocessing for co-culture scRNA-seq data

---
## Installation
Refer to [DoubletFinder GitHub](https://github.com/chris-mcginnis-ucsf/DoubletFinder) and [Symphony GitHub](https://github.com/immunogenomics/symphony).

---
## Pipeline Steps

1. **QC** (`QC_0`)
- The scripts here runs quality control on the Seurat objects created from cellranger outputs.

2. **Doublet detection** (`remove_doublets_1`)
- The scripts here runs DoubletFinder. `pK_scans_0` helps select the optimal `pK`, and `detect_doublets_1` runs DoubletFinder with the selected paramters.

3. **Seprate macrophages and fibroblasts** (`separate_fibroblasts_macrophages_2`)
- The scripts here separate macrophages and fibroblasts using canonical marker gene expression prior to reference mapping.

4. **Reference mapping with Symphony** (`reference_mapping_3`)
- The scripts here runs Symphony using the scRNA-seq atlas from Zhang et al. 2023 as the reference.
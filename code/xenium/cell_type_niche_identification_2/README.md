# Identify cell type niches in Xenium In Situ data based on the cell type identities of spatial neighbors

---

## Pipeline Steps

1. **Seurat BuildNicheAssay** (`BuildNicheAssay0`)
- The scripts here summarize the cell type profiles of the nearest 25 spatial neighbors of each cell using Seurat's BuildNicheAssay .

2. **K-means clustering to identify niches** (`KMeans1`)
- The scripts here run K-means on the niche assays from all samples to identify spatial cell type niches.
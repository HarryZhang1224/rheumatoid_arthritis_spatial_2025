# Niche-DE analysis on Xenium In Situ data

---
## Installation
Refer to [Niche-DE GitHub](https://github.com/kaishumason/NicheDE).

---
## Pipeline Steps

1. **Run NicheDE on all subtypes (combine M2 and M5)** (`M2_M5_combined_0`)
- The scripts here runs Niche-DE on all annotated subtypes to reproduce the Niche-DE related figures in figure 1. `utils/NicheDECustomGamma.R` has a custom script that allows different gamma thresholds for different cell types.
library(Seurat)
library(DoubletFinder)
library(qs)
library(glue)

sample_id <- commandArgs(trailingOnly = TRUE)
obj <- qread(glue("/data1/deyk/harry/RA_Xenium/data/Ian_old_scRNA_07062025/seurat/{sample_id}.qs"))
sweep.res <- paramSweep(obj, PCs=1:10, sct=TRUE, num.cores=10)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
res <- list(res=sweep.res,
            stats=sweep.stats,
            bcmvn=bcmvn)
qsave(res, glue("/data1/deyk/harry/RA_Xenium/results/DoubletFinder/{sample_id}.qs"))
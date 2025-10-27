#!/usr/bin/env Rscript

# -------------------------------
# Load environment and libraries
# -------------------------------
setwd("/lab-share/IM-Gutierrez-e2/Public/gvillanueva/00_CD4/analysis/")
.libPaths("/home/ch258010/miniconda3/envs/condaR441/lib/R/library")

data_lib <- c("devtools","tidyverse","msigdbr","fgsea","forcats","purrr",
              "readr","tibble","stringr","data.table")
Seurat_lib <- c("Seurat","harmony","sctransform","MAST")

print("Data wrangling packages")
lapply(data_lib, function(x) suppressMessages(suppressWarnings(require(x, character.only = TRUE))))
print("Single cell data management")
lapply(Seurat_lib, function(x) suppressMessages(suppressWarnings(require(x, character.only = TRUE))))

set.seed("42")

temp_rec <- readRDS("2025-10-24_SSc_CD4_so_process_fnl_v2.rds")

sca <- as.SingleCellExperiment(temp_rec, assay = "RNA")

filename <- paste0(format(Sys.Date(), "%Y-%m-%d_SSc_CD4_asSCexperiment_SO.rds"))

saveRDS(sca,filename)

logcounts_matrix <- as.matrix(assay(sca, "logcounts"))

# -------------------------------
# Run MAST
# -------------------------------

sca_mast <- FromMatrix(
  exprsArray = logcounts_matrix,
  cData = colData(sca),
  fData = data.frame(primerid = rownames(sca))
)

filename.2 <- paste0(format(Sys.Date(), "%Y-%m-%d_SSc_CD4_MASTmatrix_SO.rds"))

saveRDS(sca_mast,filename.2)

sca_mast<-relevel(sca_mast,"CTRL")
colData(sca_mast)$condition <- factor(colData(sca_mast)$condition, levels = c("CTRL", "SSc"))

# -------------------------------
# Run correlation
# -------------------------------
       
zlmCond <- zlm(~ condition, sca_mast)
summaryCond <- summary(zlmCond, doLRT = 'conditionSSc')
summaryDt <- summaryCond$datatable

saveRDS(zlmCond,paste0(format(Sys.Date(), "%Y-%m-%d_SSc_CD4_MAST_res.rds")))
write.table(summaryCond,"mast_pseudo_res_summary.txt")
write.table(summaryDt,"mast_pseudo_res_summary_conditions.txt")

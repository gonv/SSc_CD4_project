#!/usr/bin/env Rscript

# -------------------------------
# Load environment and libraries
# -------------------------------
.libPaths("/home/ch258010/miniconda3/envs/condaR441/lib/R/library")

data_lib <- c("devtools","tidyverse","data.table")
Seurat_lib <- c("Seurat")

# Load all packages quietly
message("Loading required packages...")
lapply(c(data_lib, Seurat_lib), function(x)
  suppressMessages(suppressWarnings(require(x, character.only = TRUE)))
)

temp_rec <- readRDS("/lab-share/IM-Gutierrez-e2/Public/gvillanueva/00_CD4/analysis/2025-10-24_SSc_CD4_so_process_fnl.rds")

temp_rec <- NormalizeData(temp_rec, normalization.method = "LogNormalize", scale.factor = 1e4)

temp_rec <- ScaleData(temp_rec,
                      vars.to.regress = c("nCount_RNA", "percent_mito", "S_score", "G2M_score"),
                      features = rownames(temp_rec))

saveRDS(temp_rec,"/lab-share/IM-Gutierrez-e2/Public/gvillanueva/00_CD4/analysis/2025-10-24_SSc_CD4_so_process_fnl_v2.rds")

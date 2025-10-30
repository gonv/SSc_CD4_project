#!/usr/bin/env Rscript

# -------------------------------
# Load environment and libraries
# -------------------------------
setwd("/lab-share/IM-Gutierrez-e2/Public/gvillanueva/00_CD4/analysis/")
.libPaths("/home/ch258010/miniconda3/envs/condaR441/lib/R/library")

Seurat_lib <- c("Seurat","SingleCellExperiment")
print("Single cell data management")
lapply(Seurat_lib, function(x) suppressMessages(suppressWarnings(require(x, character.only = TRUE))))

# -------------------------------
# Creating SCA object
# -------------------------------

# 1. Load SO object
print("Step 1")
temp_rec <- readRDS("2025-10-24_SSc_CD4_so_process_fnl_v2.rds")

# 2. "RNA" as default assay
print("Step 2")
DefaultAssay(temp_rec) <- "RNA"

# 3. Extracting matrix
print("Step 3")
counts_mat <- GetAssayData(temp_rec, assay = "RNA", layer = "counts")
logcounts_mat <- GetAssayData(temp_rec, assay = "RNA", layer = "data")
scaled_mat <- GetAssayData(temp_rec, assay = "RNA", layer = "scale.data")

# 4. Metadata
print("Step 4")
cell_metadata <- temp_rec@meta.data

# 5. Metadata by gene
print("Step 5")
gene_metadata <- data.frame(
  primerid = rownames(counts_mat),
  row.names = rownames(counts_mat)
)

# 6. Building SingleCellExperiment object
print("Step 6")
sca <- SingleCellExperiment(
  assays = list(
    counts = counts_mat,
    logcounts = logcounts_mat,
    scaledata = scaled_mat
  ),
  colData = cell_metadata,
  rowData = gene_metadata
)

# 7. Adding the reductions
print("Step 7")
if ("pca" %in% names(temp_rec@reductions)) {
  reducedDim(sca, "PCA") <- temp_rec@reductions$pca@cell.embeddings
}
if ("umap" %in% names(temp_rec@reductions)) {
  reducedDim(sca, "UMAP") <- temp_rec@reductions$umap@cell.embeddings
}

# 8. Saving
print("Step 8")
filename.1 <- paste0(format(Sys.Date(), "%Y-%m-%d_SSc_CD4_SingleCellExperiment_FULL.rds"))
saveRDS(sca, filename.1)

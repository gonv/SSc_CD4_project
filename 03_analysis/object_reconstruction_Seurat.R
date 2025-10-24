#!/usr/bin/env Rscript
# ===============================================================
# Reconstruct a Seurat object from exported Scanpy (AnnData) files
# Author: Gonzalo Villanueva-Martin
# Date: 2025-10-24
# Description:
#   This script rebuilds a Seurat object using exported data 
#   (counts, metadata, embeddings, and neighbor graphs) 
#   from a Scanpy h5ad file.
# ===============================================================

# -------------------------------
# Load environment and libraries
# -------------------------------
setwd("/lab-share/IM-Gutierrez-e2/Public/gvillanueva/00_CD4/analysis/")
.libPaths("/home/ch258010/miniconda3/envs/condaR441/lib/R/library")

# Packages
data_lib <- c("devtools","tidyverse","msigdbr","fgsea","forcats","purrr",
              "readr","tibble","stringr","data.table")
Seurat_lib <- c("Seurat","harmony","sctransform")
Single_cell <- c("Matrix","parallel","magrittr","umap","presto","ggtree","aricode","genesorteR")
ploting_lib <- c("pheatmap","RColorBrewer","ggridges","patchwork","cowplot","viridis","gridExtra","ggpubr")

# Load all packages quietly
message("Loading required packages...")
lapply(c(data_lib, Seurat_lib, Single_cell, ploting_lib), function(x)
  suppressMessages(suppressWarnings(require(x, character.only = TRUE)))
)

set.seed(42)

# -------------------------------
# Load metadata and counts
# -------------------------------
message("Loading metadata and raw counts...")

meta_ssc <- read.csv("../sc_ob/ssc_CD4_process_metadata.csv")
counts <- readMM("../sc_ob/ssc_CD4_counts_raw.mtx")

genes  <- read.csv("../sc_ob/ssc_CD4_features.tsv", sep = "\t", header = TRUE)
cells  <- read.csv("../sc_ob/ssc_CD4_barcodes.tsv", sep = "\t", header = TRUE)

# Ensure proper dimensions and names
rownames(cells) <- cells$X
cells <- cells %>% select(-X)

colnames(counts) <- genes$gene_symbol
rownames(counts) <- rownames(cells)

rownames(meta_ssc) <- meta_ssc$X
meta_ssc <- meta_ssc %>% select(-X)

# Create Seurat object
message("Creating Seurat object...")
temp_rec <- CreateSeuratObject(
  counts = t(counts),
  meta.data = meta_ssc,
  assay = "RNA",
  project = "SSc"
)

# -------------------------------
# Load neighborhood graphs
# -------------------------------
message("Loading connectivities and distances...")

conn <- readMM("../sc_ob/ssc_CD4_neighbors_connectivities.mtx")
dist <- readMM("../sc_ob/ssc_CD4_neighbors_distances.mtx")

# Store graph data into Seurat object
temp_rec@graphs$RNA_snn <- conn
temp_rec@neighbors$weighted.nn.dist <- dist

# -------------------------------
# Load PCA loadings and embeddings
# -------------------------------
message("Loading PCA loadings and embeddings...")

# PCA loadings (per gene)
pca_ssc <- read.csv("../sc_ob/ssc_CD4_pca_loadings.csv")
rownames(pca_ssc) <- pca_ssc$gene_symbol
pca_ssc <- pca_ssc %>% select(-gene_symbol)
colnames(pca_ssc) <- paste0("PC_", seq_len(ncol(pca_ssc)))

# PCA embeddings (per cell)
pca_emb <- read.csv("../sc_ob/ssc_CD4_pca_embeddings.csv")
rownames(pca_emb) <- pca_emb$X
pca_emb <- pca_emb %>% select(-X)
colnames(pca_emb) <- paste0("PC_", seq_len(ncol(pca_emb)))

pca_emb <- as.matrix(pca_emb)
storage.mode(pca_emb) <- "numeric"

# -------------------------------
# Load UMAP embeddings
# -------------------------------
message("Loading UMAP embeddings...")

umap_emb <- read.csv("../sc_ob/ssc_CD4_umap_embeddings.csv")
rownames(umap_emb) <- umap_emb$X
umap_emb <- umap_emb %>% select(-X)

umap_emb <- as.matrix(umap_emb)
storage.mode(umap_emb) <- "numeric"

# -------------------------------
# Create DimReduc objects
# -------------------------------
message("Integrating PCA and UMAP reductions...")

# PCA reduction
temp_rec[["pca"]] <- CreateDimReducObject(
  embeddings = pca_emb,
  key = "PC_",
  assay = DefaultAssay(temp_rec)
)

# Attach PCA loadings (features)
temp_rec@reductions$pca@feature.loadings <- as.matrix(pca_ssc)

# UMAP reduction
temp_rec[["umap"]] <- CreateDimReducObject(
  embeddings = umap_emb,
  key = "UMAP_",
  assay = DefaultAssay(temp_rec)
)

# -------------------------------
# Save final Seurat object
# -------------------------------
message("Saving reconstructed Seurat object...")

saveRDS(temp_rec, "2025-10-24_SSc_CD4_so_process_fnl.rds")

message("Reconstruction completed successfully.")

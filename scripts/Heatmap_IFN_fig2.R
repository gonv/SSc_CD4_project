#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyr)
  library(reshape2)
  library(tibble)
  library(openxlsx)
  library(GhibliBrewer)
  library(dplyr)
  library(purrr)
  library(pheatmap)
  library(viridis)
  library(readr)
  library(stringr)
  library(patchwork)
  library(fgsea)
  library(msigdbr)
  library(cowplot)
  library(ggh4x)
  library(Matrix)
})

source("scripts/functions.R")

# -----------------------------
# Loading input files
# -----------------------------

meta <- read.csv("data/20260124_SSc_CD4_meta_v3.csv", row.names = 1)
meta$clusters <- as.character(meta$R09)
rownames(meta) <- meta$cell
exprs_norm <- readRDS("ShinyCD4/data/20260124_SSc_CD4_counts_norm.rds")

# -----------------------------
# IFN gene sets
# -----------------------------

IFNag <- c(
  "IFITM1","TRAFD1","CD74","CASP1","PLSCR1","PSMB9","PSME1","PSME2","IFI35",
  "UBE2L6","B2M","LY6E","PSMB8","IFI30","IL7","TAP1","RTP4","IFITM2",
  "PARP9","PARP12","HELZ2","USP18","RNF31","GBP2","OAS1","BST2","LAMP3",
  "MVB12A","UBA7","TRIM21","IFI44","LPAR6","LGALS3BP","STAT2","SP110",
  "NUB1","MX1","NMI","DDX60","TDRD7","CASP8",
  "RBCK1","PTPN6","PSMB10","IFNAR2","CASP4","MVP","PIM1","ZBP1",
  "VAMP5","VAMP8","TNFSF10","OAS2","MYD88","JAK2","HLA-B","ITGB7",
  "STAT1","PSMA2","PML","FAS","STAT3","IL15RA","IL18BP","GCH1",
  "HLA-DRB1","HLA-A","SRI","PSMB2","CD274","ARID5B"
)

# -----------------------------
# Matrix by cluster – IFNa
# -----------------------------

# Keep only genes present in the expression matrix
genes_present <- IFNa_genes[IFNa_genes %in% rownames(exprs_norm)]

# Select SSc cells
cells_ssc <- meta$condition == "SSc"

# Subset expression matrix
exprs_sub <- exprs_norm[genes_present, cells_ssc, drop = FALSE]

# Subset metadata
meta_ssc <- meta[cells_ssc, ]

clusters <- factor(meta_ssc$clusters)

# Create sparse design matrix (cells x clusters)
design <- sparse.model.matrix(~ 0 + clusters)

# genes x cells  %*%  cells x clusters  → genes x clusters
exprs_cluster_sum <- exprs_sub %*% design

cluster_sizes <- colSums(design)
cluster_order <- paste0("clusters", 0:11)

# Compute mean expression per cluster
exprs_cluster_mean <- t(t(exprs_cluster_sum) / cluster_sizes)
exprs_cluster_mean <- as.matrix(exprs_cluster_mean)
exprs_cluster_mean <- exprs_cluster_mean[, cluster_order]

# Row-wise scaling (Z-score per gene)
scale_rows <- function(x) t(scale(t(x)))

exprs_scaled <- scale_rows(exprs_cluster_mean)

# Cap extreme values
exprs_scaled[exprs_scaled > 2] <- 2
exprs_scaled[exprs_scaled < -2] <- -2

exprs_scaled <- exprs_scaled[, cluster_order]

# Heatmap
pheatmap(
  t(exprs_scaled),
  color = viridis(100),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_col = 8,
  angle_col = 45,
  treeheight_row = 0,
  treeheight_col = 0,
  fontsize = 10,
  border_color = "black"
)

# -----------------------------
# IFNg by cluster
# -----------------------------

# Keep only genes present
genes_present <- IFNg_genes[IFNg_genes %in% rownames(exprs_norm)]

# Select SSc cells
cells_ssc <- meta$condition == "SSc"

# Subset expression matrix
exprs_sub <- exprs_norm[genes_present, cells_ssc, drop = FALSE]

# Subset metadata
meta_ssc <- meta[cells_ssc, ]

clusters <- factor(meta_ssc$clusters)

# Create sparse design matrix
design <- sparse.model.matrix(~ 0 + clusters)

# Compute cluster-level means
exprs_cluster_sum <- exprs_sub %*% design
cluster_sizes <- colSums(design)
cluster_order <- paste0("clusters", 0:11)

exprs_cluster_mean <- t(t(exprs_cluster_sum) / cluster_sizes)
exprs_cluster_mean <- as.matrix(exprs_cluster_mean)
exprs_cluster_mean <- exprs_cluster_mean[, cluster_order]

# Scale per gene
exprs_scaled <- scale_rows(exprs_cluster_mean)

# Cap extreme values
exprs_scaled[exprs_scaled > 2] <- 2
exprs_scaled[exprs_scaled < -2] <- -2

exprs_scaled <- exprs_scaled[, cluster_order]

# Heatmap
pheatmap(
  t(exprs_scaled),
  color = viridis(100),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_col = 8,
  angle_col = 45,
  treeheight_row = 0,
  treeheight_col = 0,
  fontsize = 10,
  border_color = "black"
)

# -----------------------------
# Expression by cluster + condition
# -----------------------------

# Keep only genes present
genes_present <- IFNag[IFNag %in% rownames(exprs_norm)]
exprs_sub <- exprs_norm[genes_present, , drop = FALSE]

# Define desired ordering
cluster_order <- 0:11
condition_order <- c("CTRL", "SSc")

desired_levels <- with(
  expand.grid(condition = condition_order,
              cluster = cluster_order),
  paste0("cluster", cluster, "_", condition)
)

# Create combined cluster-condition factor
meta$cluster_condition <- paste0(
  "cluster", meta$clusters, "_", meta$condition
)

group_factor <- factor(meta$cluster_condition,
                       levels = desired_levels)

# Sparse design matrix (cells x groups)
design <- sparse.model.matrix(~ 0 + group_factor)

# Compute mean per group
exprs_sum <- exprs_sub %*% design
group_sizes <- colSums(design)
exprs_mean <- t(t(exprs_sum) / group_sizes)

exprs_mean <- as.matrix(exprs_mean)
colnames(exprs_mean) <- gsub("group_factor", "", colnames(exprs_mean))

# Scale per gene
exprs_scaled <- scale_rows(exprs_mean)

# Cap extreme values
exprs_scaled[exprs_scaled > 2] <- 2
exprs_scaled[exprs_scaled < -2] <- -2

# Heatmap
pheatmap(
  exprs_scaled,
  color = viridis(100),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_col = 10,
  angle_col = 45,
  treeheight_row = 0,
  treeheight_col = 0,
  fontsize = 12,
  border_color = "white"
)

# -----------------------------
# Signature score heatmap
# -----------------------------

# Keep only genes present
IFNa_present <- IFNa_genes[IFNa_genes %in% rownames(exprs_norm)]
IFNg_present <- IFNg_genes[IFNg_genes %in% rownames(exprs_norm)]

# Compute per-cell scores
meta$IFNa_score <- Matrix::colMeans(exprs_norm[IFNa_present, , drop = FALSE])
meta$IFNg_score <- Matrix::colMeans(exprs_norm[IFNg_present, , drop = FALSE])

# Compute mean score per cluster and condition
score_df <- meta %>%
  group_by(clusters, condition) %>%
  summarise(
    IFNa = mean(IFNa_score),
    IFNg = mean(IFNg_score),
    .groups = "drop"
  )

# Convert to wide format
score_wide <- score_df %>%
  pivot_wider(
    names_from = condition,
    values_from = c(IFNa, IFNg),
    names_sep = "_"
  )

score_wide$clusters <- factor(score_wide$clusters,
                              levels = as.character(0:11))

score_wide <- score_wide %>%
  arrange(clusters)

score_mat <- as.matrix(score_wide[, -1])
rownames(score_mat) <- paste0("cluster_", score_wide$clusters)

score_mat <- score_mat[, c(
  "IFNa_CTRL","IFNa_SSc",
  "IFNg_CTRL","IFNg_SSc"
)]

# Column annotation
annotation_col <- data.frame(
  Signature = c("IFNa","IFNa","IFNg","IFNg"),
  Condition = c("CTRL","SSc","CTRL","SSc")
)

rownames(annotation_col) <- colnames(score_mat)

annotation_colors <- list(
  Signature = c(IFNa = "#5E548E", IFNg = "#1B998B"),
  Condition = c(CTRL = "#D8BFD8", SSc = "#C3E5C3")
)

# Final heatmap
pheatmap(
  t(score_mat),
  color = viridis(100),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  border_color = "black",
  show_rownames = FALSE,
  annotation_row = annotation_col,
  annotation_colors = annotation_colors,
  gaps_row = 2,
  angle_col = 45,
  fontsize_col = 9,
  legend = TRUE,
  annotation_legend = TRUE,
  fontsize = 10,
  annotation_names_row = FALSE,
  annotation_names_col = FALSE
)

# -----------------------------
# Th heatmap
# -----------------------------

# Select Tcm cells
cells_tcm <- meta$cell_type == "Tcm"

meta_th <- meta[cells_tcm, ]
exprs_tcm <- exprs_norm[, cells_tcm, drop = FALSE]

# Keep only genes present
genes_present <- genes_th[genes_th %in% rownames(exprs_tcm)]
exprs_sub <- exprs_tcm[genes_present, , drop = FALSE]

# Create Th + condition grouping
meta_th$Th_condition <- paste(meta_th$Th, meta_th$condition, sep = "_")

groups <- factor(meta_th$Th_condition)

design <- sparse.model.matrix(~ 0 + groups)
colnames(design) <- levels(groups)

# Compute mean per group
exprs_sum <- exprs_sub %*% design
group_sizes <- colSums(design)

valid_groups <- group_sizes > 0
exprs_sum <- exprs_sum[, valid_groups, drop = FALSE]
group_sizes <- group_sizes[valid_groups]

exprs_mean <- t(t(exprs_sum) / group_sizes)
exprs_mean <- as.matrix(exprs_mean)

# Alternating CTRL/SSc order within each Th
Th_levels <- sort(unique(meta_th$Th))
condition_order <- c("CTRL", "SSc")

desired_order <- as.vector(
  t(outer(Th_levels, condition_order,
          function(th, cond) paste0(th, "_", cond)))
)

desired_order <- intersect(desired_order, colnames(exprs_mean))
exprs_mean <- exprs_mean[, desired_order, drop = FALSE]

# Robust row scaling
scale_rows_safe <- function(mat) {
  row_means <- rowMeans(mat)
  row_sds <- apply(mat, 1, sd)
  row_sds[row_sds == 0] <- 1
  scaled <- sweep(mat, 1, row_means, "-")
  scaled <- sweep(scaled, 1, row_sds, "/")
  return(scaled)
}

exprs_scaled <- scale_rows_safe(exprs_mean)

# Cap extreme values
exprs_scaled[exprs_scaled > 2] <- 2
exprs_scaled[exprs_scaled < -2] <- -2

# Heatmap
pheatmap(
  exprs_scaled,
  color = viridis(100),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  border_color = "black",
  angle_col = 90,
  fontsize_col = 10
)
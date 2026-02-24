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

# Loading files----

meta <- read.csv("data/20260124_SSc_CD4_meta_v3.csv", row.names = 1)
meta$clusters <- as.character(meta$R09)
rownames(meta) <- meta$cell
exprs_norm <- readRDS("ShinyCD4/data/20260124_SSc_CD4_counts_norm.rds")

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

IFNa_genes <- c(
  "IFITM1","TRAFD1","CD74","CASP1","PLSCR1","PSMB9","PSME1","PSME2",
  "IFI35","UBE2L6","B2M","LY6E","PSMB8","IFI30","IL7","TAP1","RTP4",
  "IFITM2","PARP9","PARP12","HELZ2","USP18","RNF31","GBP2","OAS1",
  "BST2","LAMP3","MVB12A","UBA7","TRIM21","IFI44","LPAR6",
  "LGALS3BP","STAT2","SP110","NUB1","MX1","NMI","DDX60",
  "TDRD7","CASP8"
)

IFNg_genes <- c(
  "RBCK1","TRAFD1","PTPN6","CD74","PSMB10","CASP1","PLSCR1","PSMB9",
  "PSME1","PSME2","IFNAR2","IFI35","UBE2L6","CASP4","MVP","B2M",
  "LY6E","PSMB8","IFI30","PIM1","IL7","ZBP1","TAP1","RTP4",
  "VAMP5","VAMP8","IFITM2","PARP12","TNFSF10","OAS2","HELZ2",
  "USP18","RNF31","MYD88","BST2","JAK2","HLA-B","ITGB7",
  "STAT1","TRIM21","PSMA2","PML","FAS","IFI44","LGALS3BP",
  "STAT3","IL15RA","STAT2","IL18BP","SP110","GCH1","MX1",
  "HLA-DRB1","HLA-A","NMI","DDX60","SRI","TDRD7",
  "PSMB2","CD274","ARID5B","CASP8"
)


# matrix by cluster ----

genes_present <- IFNa_genes[IFNa_genes %in% rownames(exprs_norm)]

#exprs_sub <- exprs_norm[genes_present, , drop = FALSE]
# Índices de células SSc
cells_ssc <- meta$condition == "SSc"

# Filtrar expresión
exprs_sub <- exprs_norm[genes_present, cells_ssc, drop = FALSE]

# Filtrar metadata
meta_ssc <- meta[cells_ssc, ]

clusters <- factor(meta_ssc$clusters)
#clusters <- factor(meta$clusters)

# Crear matriz de diseño cells x clusters
design <- sparse.model.matrix(~ 0 + clusters)

# genes x cells  %*%  cells x clusters  → genes x clusters
exprs_cluster_sum <- exprs_sub %*% design

cluster_sizes <- colSums(design)
cluster_order <- paste0("clusters", 0:11)


exprs_cluster_mean <- t(t(exprs_cluster_sum) / cluster_sizes)

exprs_cluster_mean <- as.matrix(exprs_cluster_mean)
exprs_cluster_mean <- exprs_cluster_mean[, cluster_order]

scale_rows <- function(x) t(scale(t(x)))

exprs_scaled <- scale_rows(exprs_cluster_mean)
exprs_scaled[exprs_scaled > 2] <- 2
exprs_scaled[exprs_scaled < -2] <- -2

exprs_scaled <- exprs_scaled[, cluster_order]

pheatmap(
  t(exprs_scaled),
  color = viridis(100),
  
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_col = 8,
  angle_col = 45,
  
  # estética
  treeheight_row = 0,   # equivalente a show_row_dend = FALSE
  treeheight_col = 0,   # equivalente a show_column_dend = FALSE
  
  fontsize = 10,
  border_color = "black"
  
)

# IFNg by cluster----
genes_present <- IFNg_genes[IFNg_genes %in% rownames(exprs_norm)]

#exprs_sub <- exprs_norm[genes_present, , drop = FALSE]

#clusters <- factor(meta$clusters)

cells_ssc <- meta$condition == "SSc"

# Filtrar expresión
exprs_sub <- exprs_norm[genes_present, cells_ssc, drop = FALSE]

# Filtrar metadata
meta_ssc <- meta[cells_ssc, ]

clusters <- factor(meta_ssc$clusters)


# Crear matriz de diseño cells x clusters
design <- sparse.model.matrix(~ 0 + clusters)

# genes x cells  %*%  cells x clusters  → genes x clusters
exprs_cluster_sum <- exprs_sub %*% design

cluster_sizes <- colSums(design)
cluster_order <- paste0("clusters", 0:11)


exprs_cluster_mean <- t(t(exprs_cluster_sum) / cluster_sizes)

exprs_cluster_mean <- as.matrix(exprs_cluster_mean)
exprs_cluster_mean <- exprs_cluster_mean[, cluster_order]

scale_rows <- function(x) t(scale(t(x)))

exprs_scaled <- scale_rows(exprs_cluster_mean)
exprs_scaled[exprs_scaled > 2] <- 2
exprs_scaled[exprs_scaled < -2] <- -2

exprs_scaled <- exprs_scaled[, cluster_order]

pheatmap(
  t(exprs_scaled),
  color = viridis(100),
  
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_col = 8,
  angle_col = 45,
  
  # estética
  treeheight_row = 0,   # equivalente a show_row_dend = FALSE
  treeheight_col = 0,   # equivalente a show_column_dend = FALSE
  
  fontsize = 10,
  border_color = "black"
  
)

# by condition ----

### 1️⃣ Genes presentes
genes_present <- IFNag[IFNag %in% rownames(exprs_norm)]
exprs_sub <- exprs_norm[genes_present, , drop = FALSE]

### 2️⃣ Crear factor combinado cluster + condition
### Definir orden deseado
cluster_order <- 0:11
condition_order <- c("CTRL", "SSc")

desired_levels <- with(
  expand.grid(condition = condition_order,
              cluster = cluster_order),
  paste0("cluster", cluster, "_", condition)
)

### Crear factor con levels en orden correcto
meta$cluster_condition <- paste0(
  "cluster", meta$clusters, "_", meta$condition
)

group_factor <- factor(meta$cluster_condition,
                       levels = desired_levels)

### 3️⃣ Crear matriz de diseño sparse (cells x groups)
design <- sparse.model.matrix(~ 0 + group_factor)

### 4️⃣ Sumar expresión por grupo (genes x groups)
exprs_sum <- exprs_sub %*% design

### 5️⃣ Dividir por número de células por grupo
group_sizes <- colSums(design)
exprs_mean <- t(t(exprs_sum) / group_sizes)

### 6️⃣ Convertir a matriz densa (ahora es pequeña)
exprs_mean <- as.matrix(exprs_mean)

exprs_mean <- exprs_mean[, desired_order]

### 8️⃣ Limpiar nombres columnas
colnames(exprs_mean) <- gsub("group_factor", "", colnames(exprs_mean))

### 9️⃣ Escalar por gen (Z-score)
scale_rows <- function(x) t(scale(t(x)))
exprs_scaled <- scale_rows(exprs_mean)

### 🔟 Capar extremos
exprs_scaled[exprs_scaled > 2] <- 2
exprs_scaled[exprs_scaled < -2] <- -2

### 1️⃣1️⃣ Heatmap
pheatmap(
  exprs_scaled,
  color = viridis(100),
  
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_col = 10,
  angle_col = 45,
  
  # estética
  treeheight_row = 0,   # equivalente a show_row_dend = FALSE
  treeheight_col = 0,   # equivalente a show_column_dend = FALSE
  
  fontsize = 12,
  border_color = "white"
)

annotation_col <- data.frame(
  condition = ifelse(grepl("SSc", colnames(exprs_scaled)), "SSc", "CTRL")
)
rownames(annotation_col) <- colnames(exprs_scaled)

pheatmap(
  exprs_scaled,
  color = viridis(100),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_col = 10,
  angle_col = 45,
  
  # estética
  treeheight_row = 0,   # equivalente a show_row_dend = FALSE
  treeheight_col = 0,   # equivalente a show_column_dend = FALSE
  
  fontsize = 12,
  border_color = "white",
  annotation_col = annotation_col
)

# Heatmap by cluster----
### 1️⃣ Genes presentes
IFNa_present <- IFNa_genes[IFNa_genes %in% rownames(exprs_norm)]
IFNg_present <- IFNg_genes[IFNg_genes %in% rownames(exprs_norm)]

### 2️⃣ Score por célula
meta$IFNa_score <- Matrix::colMeans(exprs_norm[IFNa_present, , drop = FALSE])
meta$IFNg_score <- Matrix::colMeans(exprs_norm[IFNg_present, , drop = FALSE])

### 3️⃣ Media por cluster + condición
score_df <- meta %>%
  group_by(clusters, condition) %>%
  summarise(
    IFNa = mean(IFNa_score),
    IFNg = mean(IFNg_score),
    .groups = "drop"
  )

### 4️⃣ Formato ancho
score_wide <- score_df %>%
  pivot_wider(
    names_from = condition,
    values_from = c(IFNa, IFNg),
    names_sep = "_"
  )

### 5️⃣ Ordenar clusters y renombrar
score_wide$clusters <- factor(score_wide$clusters,
                              levels = as.character(0:11))

score_wide <- score_wide %>%
  arrange(clusters)

score_mat <- as.matrix(score_wide[, -1])
rownames(score_mat) <- paste0("cluster_", score_wide$clusters)

### 6️⃣ Orden columnas explícito
score_mat <- score_mat[, c(
  "IFNa_CTRL","IFNa_SSc",
  "IFNg_CTRL","IFNg_SSc"
)]

### 7️⃣ Crear anotación lateral (para columnas)
annotation_col <- data.frame(
  Signature = c("IFNa","IFNa","IFNg","IFNg"),
  Condition = c("CTRL","SSc","CTRL","SSc")
)

rownames(annotation_col) <- colnames(score_mat)

annotation_colors <- list(
  Signature = c(IFNa = "#5E548E", IFNg = "#1B998B"),
  Condition = c(CTRL = "#D8BFD8", SSc = "#C3E5C3")
)

score_mat_plot <- score_mat
colnames(score_mat_plot) <- NULL

### 2️⃣ Separación visual entre IFNa e IFNg
# (tenemos 4 filas: IFNa_CTRL, IFNa_SSc, IFNg_CTRL, IFNg_SSc)
gaps_row <- 2   # separa después de IFNa

### 3️⃣ Heatmap final
pheatmap(
  t(score_mat_plot),
  color = viridis(100),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  border_color = "black",
  show_rownames = FALSE,
  
  annotation_row = annotation_col,
  annotation_colors = annotation_colors,
  
  gaps_row = gaps_row,
  
  angle_col = 45,          # ligera rotación
  fontsize_col = 9,
  
  legend = TRUE,
  legend_breaks = pretty(range(score_mat_plot), n = 5),
  
  annotation_legend = TRUE,
  
  legend_labels = NULL,
  
  fontsize = 10,
  
  annotation_names_row = FALSE,
  
  annotation_names_col = FALSE
)


# Heatmap Th----
### 1️⃣ Filtrar Tcm
cells_tcm <- meta$cell_type == "Tcm"

meta_th <- meta[cells_tcm, ]
exprs_tcm <- exprs_norm[, cells_tcm, drop = FALSE]

### 2️⃣ Genes presentes
genes_present <- genes_th[genes_th %in% rownames(exprs_tcm)]
exprs_sub <- exprs_tcm[genes_present, , drop = FALSE]

### 3️⃣ Crear combinación Th + condition
meta_th$Th_condition <- paste(meta_th$Th, meta_th$condition, sep = "_")

### 4️⃣ Crear diseño sparse SIN prefijo
groups <- factor(meta_th$Th_condition)

design <- sparse.model.matrix(~ 0 + groups)
colnames(design) <- levels(groups)  # eliminar prefijo automático

### 5️⃣ Media por grupo
exprs_sum <- exprs_sub %*% design
group_sizes <- colSums(design)

# eliminar posibles grupos vacíos
valid_groups <- group_sizes > 0
exprs_sum <- exprs_sum[, valid_groups, drop = FALSE]
group_sizes <- group_sizes[valid_groups]

exprs_mean <- t(t(exprs_sum) / group_sizes)
exprs_mean <- as.matrix(exprs_mean)

### 6️⃣ Orden alternante CTRL / SSc dentro de cada Th
Th_levels <- sort(unique(meta_th$Th))
condition_order <- c("CTRL", "SSc")

desired_order <- as.vector(
  t(outer(Th_levels, condition_order,
          function(th, cond) paste0(th, "_", cond)))
)

desired_order <- intersect(desired_order, colnames(exprs_mean))

exprs_mean <- exprs_mean[, desired_order, drop = FALSE]

### 7️⃣ Escalado robusto por gen
scale_rows_safe <- function(mat) {
  row_means <- rowMeans(mat)
  row_sds <- apply(mat, 1, sd)
  row_sds[row_sds == 0] <- 1
  
  scaled <- sweep(mat, 1, row_means, "-")
  scaled <- sweep(scaled, 1, row_sds, "/")
  return(scaled)
}

exprs_scaled <- scale_rows_safe(exprs_mean)

exprs_scaled[exprs_scaled > 2] <- 2
exprs_scaled[exprs_scaled < -2] <- -2

### 8️⃣ Heatmap
pheatmap(
  exprs_scaled,
  color = viridis(100),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  border_color = "black",
  angle_col = 90,
  fontsize_col = 10
)


#

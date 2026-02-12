#!/usr/bin/env Rscript

# sc info ####
#dea_pb <- read.csv("data/dea_pseudobulk_condition_v5.csv")
exprs_norm <- readRDS("data/20260124_SSc_CD4_counts_norm.rds")
#saveRDS(exprs_norm, file = "data/exprs_norm.rds", compress = "xz")
#counts <- readRDS("data/exprs_norm.rds")

meta <- read.csv("data/20260124_SSc_CD4_meta_v2.csv")
meta$clusters <- as.character(meta$R09)
meta$cell <- as.character(meta$X)
IR_V.cols <- colnames(meta[,grep("IR_V*", colnames(meta))])
ir_vdj <- meta %>% 
  select("X","Sample","condition","cell_type",IR_V.cols)
meta <- meta %>% 
  select(!IR_V.cols, )

save(list = c("meta","ir_vdj"),file = "SScCD4_data.RData")

files_cluster <- list.files("data/dea_pseudobulk/cluster/")

# DEA ####

clean_dea_name <- function(filename, prefix) {
  basename <- tools::file_path_sans_ext(basename(filename))
  cleaned <- sub("dea_pseudobulk_", "", basename)
}

clean_dea_CT_name <- function(filename, prefix) {
  basename <- tools::file_path_sans_ext(basename(filename))
  cleaned <- sub("dea_pseudobulk_CT_", "", basename)
}

# 1. Ruta principal
base_path <- "data/dea_pseudobulk"

# 2. Leer archivo principal
dea_main <- read.csv(file.path(base_path, "dea_pseudobulk_condition_v5.csv"))

# 3. Leer subcarpeta de cluster
cluster_files <- list.files(file.path(base_path, "cluster"), pattern = "*.csv", full.names = TRUE)
cluster_list <- lapply(cluster_files, read.csv)
names(cluster_list) <- vapply(cluster_files, clean_dea_name, character(1))

# 4. Leer subcarpeta de cell type
celltype_files <- list.files(file.path(base_path, "cell_type"), pattern = "*.csv", full.names = TRUE)
celltype_list <- lapply(celltype_files, read.csv)
names(celltype_list) <- vapply(celltype_files, clean_dea_CT_name, character(1))

# 5. Combinar todo
dea_list_ps <- c(
  Global = list(dea_main),
  cluster_list,
  celltype_list
)

# 6. Guardar
save(dea_list_ps, file = "data/dea_list_ps.RData")

# GSEA ####

# 1. Ruta principal
base_path <- "data/GSEA_TFactivity"

# 1. Función auxiliar para limpiar nombres
clean_gsea_name <- function(filename, prefix) {
  basename <- tools::file_path_sans_ext(basename(filename))
  cleaned <- sub("^2026-01-\\d{2}_SSc_CD4_GSEA_", "", basename)
}

clean_Reactome_name <- function(filename, prefix) {
  basename <- tools::file_path_sans_ext(basename(filename))
  cleaned <- sub("^2026-01-\\d{2}_SSc_CD4_GSEA_", "", basename)
  cleaned <- sub("C2", "Reactome", cleaned)
}

# 2. Crear las listas con nombres limpios

H_files <- list.files(file.path(base_path, "Hallmark"), pattern = "*.csv", full.names = TRUE)
H_list <- lapply(H_files, read.csv)
gseaHallmark_main <- read.csv("data/GSEA_TFactivity/2026-01-31_SSc_CD4_GSEA_pseudobulk_Hallmark.csv")
names(H_list) <- vapply(H_files, clean_gsea_name, character(1))
H_list[["Global_Hallmark"]] <- gseaHallmark_main

C2_files <- list.files(file.path(base_path, "C2"), pattern = "*.csv", full.names = TRUE)
C2_list <- lapply(C2_files, read.csv)
gseaC2_main <- read.csv("data/GSEA_TFactivity/2026-01-31_SSc_CD4_GSEA_pseudobulk_C2.csv")
names(C2_list) <- vapply(C2_files, clean_Reactome_name, character(1))
C2_list[["Global_Reactome"]] <- gseaC2_main

TF_files <- list.files(file.path(base_path, "TFactivity"), pattern = "*.csv", full.names = TRUE)
TF_list <- lapply(TF_files, read.csv)
TF_main <- read.csv("data/GSEA_TFactivity/2026-01-31_SSc_CD4_GSEA_pseudobulk_TFact.csv")
names(TF_list) <- vapply(TF_files, clean_gsea_name, character(1))
TF_list[["Global_TFactivity"]] <- TF_main

# 3. Combinar todas las listas
gseaTF_list <- c(H_list, C2_list, TF_list)


# 6. Guardar
save(gseaTF_list, file = "data/gseaTF_list.RData")

# DEA version 2####

clean_dea_name <- function(filename, prefix) {
  basename <- tools::file_path_sans_ext(basename(filename))
  cleaned <- sub("scRNAseq_SScCD4_", "", basename)
  cleaned <- sub("_DEG_R09_byCON_1500HVG_20PC_it92", "", cleaned)
}

clean_dea_CT_name <- function(filename, prefix) {
  basename <- tools::file_path_sans_ext(basename(filename))
  cleaned <- sub("scRNAseq_SScCD4_", "", basename)
  cleaned <- sub("_DEG_R09_byCON_1500HVG_20PC_it92", "", cleaned)
}

scanpy_to_dea <- function(df) {
  out <- data.frame(
    gene   = df$SSc_names,
    logFC  = df$SSc_logfoldchanges,
    pval   = df$SSc_pvals,
    padj   = df$SSc_pvals_adj,
    zscore   = df$SSc_scores,
    #cluster = cluster_id,
    stringsAsFactors = FALSE
  )
  
  out <- out[!is.na(out$gene), ]
  return(out)
}

# 1. Ruta principal
base_path <- "data/dea_sc/"

# 2. Leer archivo principal
dea_main <- read.csv(file.path(base_path, "../dea_pseudobulk/dea_pseudobulk_condition_v5.csv"))

# 3. Leer subcarpeta de cluster
cluster_files <- list.files(file.path(base_path, "cluster"), pattern = "*.csv", full.names = TRUE)
cluster_list <- lapply(cluster_files, read.csv)
names(cluster_list) <- vapply(cluster_files, clean_dea_name, character(1))

cluster_id <- names(cluster_list)

dea_clusters <- lapply(names(cluster_list), function(cl) {
  scanpy_to_dea(cluster_list[[cl]])
})

names(dea_clusters) <- names(cluster_list)

# 4. Leer subcarpeta de cell type
celltype_files <- list.files(file.path(base_path, "cell_type"), pattern = "*.csv", full.names = TRUE)
celltype_list <- lapply(celltype_files, read.csv)
names(celltype_list) <- vapply(celltype_files, clean_dea_CT_name, character(1))

dea_ct <- lapply(names(celltype_list), function(cl) {
  scanpy_to_dea(celltype_list[[cl]])
})

names(dea_ct) <- names(celltype_list)

# 5. Combinar todo
dea_list <- c(
  Global = list(dea_main),
  dea_clusters,
  dea_ct
)

# 6. Guardar
save(dea_list, file = "data/dea_list_2.RData")

save(list = c("meta","dea_list","dea_list_ps","gseaTF_list","exprs_norm"),file = "data/SScCD4_data_v2.RData")

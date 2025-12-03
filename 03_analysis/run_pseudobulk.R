#!/usr/bin/env Rscript

## ============================================================
## Gonzalo Villanueva-Martin
## Pseudobulk differential expression (CTRL vs SSc) with LM
## Adjusted for donor and library size (cell_umi)
## ============================================================

setwd("/lab-share/IM-Gutierrez-e2/Public/gvillanueva/00_CD4/analysis/")

.libPaths("/home/ch258010/miniconda3/envs/condaR441/lib/R/library")
set.seed(42)

## ------------------------------------------------------------
## Load libraries
## ------------------------------------------------------------
data_lib <- c("devtools","tidyverse","forcats","purrr","readr",
              "tibble","stringr","data.table","DT")
pathway_lib <- c("msigdbr","fgsea")
RNAseq_lib <- c("edgeR","lme4")
Seurat_lib <- c("Seurat","harmony","sctransform","MAST","scater")
Single_cell <- c("Matrix","parallel","magrittr","umap","presto",
                 "ggtree","aricode","genesorteR")
ploting_lib <- c("pheatmap","RColorBrewer","ggridges","patchwork",
                 "cowplot","viridis","gridExtra","ggpubr","NMF",
                 "purrr","ggrepel")

message("Data wrangling packages")
lapply(data_lib, function(x) suppressMessages(suppressWarnings(require(x, character.only = TRUE))))

message("Pathway analysis packages")
lapply(pathway_lib, function(x) suppressMessages(suppressWarnings(require(x, character.only = TRUE))))      

message("RNA-seq data management")
lapply(RNAseq_lib, function(x) suppressMessages(suppressWarnings(require(x, character.only = TRUE))))      

message("Single cell data management")
lapply(Seurat_lib, function(x) suppressMessages(suppressWarnings(require(x, character.only = TRUE))))
lapply(Single_cell, function(x) suppressMessages(suppressWarnings(require(x, character.only = TRUE))))     

message("Plotting data packages")
lapply(ploting_lib, function(x) suppressMessages(suppressWarnings(require(x, character.only = TRUE))))

## ------------------------------------------------------------
## Load helper functions
## ------------------------------------------------------------
message("Loading custom functions")
source("/lab-share/IM-Gutierrez-e2/Public/gvillanueva/00_CD4/scripts/someFunctions_SR.R")
source("/lab-share/IM-Gutierrez-e2/Public/Qian/B_Atlas/analysis_scripts/Some_functions.r")

## ------------------------------------------------------------
## Differential expression function (per gene)
## ------------------------------------------------------------
diffexp_pseudobulk <- function(gexp, gname, group, donor, cell_umi) {
  # gexp: numeric vector of log2(CPM+1) for one gene across samples
  # group: factor with levels c("CTRL","SSc")
  # donor: factor indicating donor ID per sample
  # cell_umi: numeric covariate (log2 total UMI per sample)
  
  df <- data.frame(
    gexp     = as.numeric(gexp),
    group    = group,
    donor    = donor,
    cell_umi = cell_umi
  )
  
  # Full model: condition + donor + library size
  m_full <- lm(gexp ~ group + donor + cell_umi, data = df)
  
  # Null model: donor + library size only
  m_null <- lm(gexp ~ donor + cell_umi, data = df)
  
  # Likelihood ratio test between null and full model
  m_lrt <- anova(m_null, m_full, test = "LRT")
  
  # Extract coefficient for the disease effect (SSc vs CTRL)
  coef_summary <- summary(m_full)$coefficients
  if (!"groupSSc" %in% rownames(coef_summary)) {
    # In case the model drops the contrast (rare, but just in case)
    return(data.frame(
      Gene       = gname,
      Estimate   = NA,
      StdError   = NA,
      tValue     = NA,
      Pvalue     = NA,
      LRT_ChiSq  = NA,
      LRT_Pvalue = NA
    ))
  }
  
  beta <- coef_summary["groupSSc", ]
  
  data.frame(
    Gene       = gname,
    Estimate   = beta["Estimate"],
    StdError   = beta["Std. Error"],
    tValue     = beta["t value"],
    Pvalue     = beta["Pr(>|t|)"],
    LRT_ChiSq  = m_lrt$`Sum of Sq`[2],
    LRT_Pvalue = m_lrt$`Pr(>Chi)`[2]
  )
}

## ------------------------------------------------------------
## Load Seurat object and prepare counts / metadata
## ------------------------------------------------------------
message("Loading Seurat object")
temp_rec <- readRDS("2025-10-24_SSc_CD4_so_process_fnl_v3.rds")

# Raw counts (genes x cells)
counts <- temp_rec@assays$RNA@layers$counts

# Metadata
meta_ssc <- temp_rec@meta.data
meta_ssc$cell <- rownames(meta_ssc)

## ------------------------------------------------------------
## Generate pseudobulk matrix
## ------------------------------------------------------------
message("Collapsing counts into pseudobulk samples")

# Ensure column order matches metadata
stopifnot(all(colnames(counts) == meta_ssc$cell))

all_collapse <- collapse_counts(
  counts,
  meta_ssc,
  c("R09", "Sample", "condition")
)
# all_collapse$counts_mat: genes x pseudobulk-samples
# all_collapse$meta_data: data.frame with R09, Sample, condition, sample_id

## ------------------------------------------------------------
## Gene filtering: remove mitochondrial / ribosomal / MALAT1 / MIR-
## ------------------------------------------------------------
genes_exclude <- grep("^MT-|^RPL|^RPS|MALAT1|MIR-", 
                      rownames(all_collapse$counts_mat), 
                      value = TRUE)
keep_genes <- setdiff(rownames(all_collapse$counts_mat), genes_exclude)

counts_filtered <- all_collapse$counts_mat[keep_genes, ]

## ------------------------------------------------------------
## Filter lowly expressed genes with edgeR::filterByExpr
## ------------------------------------------------------------
message("Filtering lowly expressed genes using edgeR::filterByExpr")

y <- DGEList(counts = counts_filtered)
group <- factor(all_collapse$meta_data$condition, levels = c("CTRL", "SSc"))
design <- model.matrix(~ group)

keep <- filterByExpr(y, design)
y <- y[keep, , keep.lib.sizes = FALSE]

d <- y$counts  # filtered count matrix (genes x pseudobulk-samples)

## ------------------------------------------------------------
## Compute log2(CPM+1) matrix for modeling
## ------------------------------------------------------------
message("Computing log2(CPM + 1) expression matrix")

tmp <- log2(edgeR::cpm(d) + 1)
tmp <- as(tmp, "dgCMatrix")  # for consistency with sparse operations
gene_ids <- rownames(tmp)

## ------------------------------------------------------------
## Build design variables for the linear model
## ------------------------------------------------------------
donor <- factor(all_collapse$meta_data$Sample)  # donor/sample ID per pseudobulk
cell_umi <- apply(all_collapse$counts_mat, 2, function(x) log2(sum(x) + 1))

# Safety checks
stopifnot(length(group)    == ncol(tmp))
stopifnot(length(donor)    == ncol(tmp))
stopifnot(length(cell_umi) == ncol(tmp))

## ------------------------------------------------------------
## Run gene-wise linear models
## ------------------------------------------------------------
message("Running gene-wise linear models (pseudobulk DE)")

results_list <- lapply(seq_len(nrow(tmp)), function(i) {
  diffexp_pseudobulk(
    gexp  = tmp[i, ],
    gname = gene_ids[i],
    group = group,
    donor = donor,
    cell_umi = cell_umi
  )
})

results <- do.call(rbind, results_list)

## ------------------------------------------------------------
## Prepare result table with FDR and labels
## ------------------------------------------------------------
results_temp <- results[, c("Gene","Estimate","StdError","tValue","Pvalue")]
colnames(results_temp) <- c("gene","Estimate", "Std_error", "t_value", "p_value")
rownames(results_temp) <- seq_len(nrow(results_temp))

results_temp <- as.data.frame(results_temp)
results_temp$log_p_value <- -log10(results_temp$p_value)
results_temp$adj.P.Val   <- p.adjust(results_temp$p_value, method = "fdr")

# Label significant results
results_temp$significance <- dplyr::case_when(
  results_temp$adj.P.Val < 0.05 & results_temp$Estimate >  0.25 ~ "Upregulated",
  results_temp$adj.P.Val < 0.05 & results_temp$Estimate < -0.25 ~ "Downregulated",
  TRUE ~ "Not_significant"
)

results_temp$delabel <- NA
results_temp$delabel[results_temp$significance != "Not_significant"] <-
  results_temp$gene[results_temp$significance != "Not_significant"]

## ------------------------------------------------------------
## Save results
## ------------------------------------------------------------
message("Saving results to CSV")
write.csv(results_temp, "pseudobulk_condition_v4.csv", row.names = FALSE)

message("Pseudobulk DE analysis completed successfully.")

#!/usr/bin/env Rscript

# ------------------------------------------------------------
# Figure 1 generation script
# CD4 T cell single-cell RNA-seq analysis in SSc
# ------------------------------------------------------------

libs <- c(
  "tidyverse",
  "data.table",
  "ggplot2",
  "ggrepel",
  "pheatmap",
  "patchwork",
  "RColorBrewer",
  "ggpubr",
  "decoupleR"
)

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(ggplot2)
  library(pheatmap)
  library(pheatmap)
  library(patchwork)
  library(RColorBrewer)
  library(ggpubr)
  library(ggrepel)
  library(decoupleR)
})

lapply(libs, function(x) suppressMessages(suppressWarnings(require(x, character.only = TRUE))))

#############################################################
# Load metadata, counts and DEG results
#############################################################

res2_ann <- read.csv("ShinyCD4/data/dea_pseudobulk_condition_v5.csv")
exprs_norm <- readRDS("ShinyCD4/data/20260124_SSc_CD4_counts_norm.rds")
meta <- read.csv("data/20260124_SSc_CD4_meta_v3.csv", row.names = 1)

#############################################################
# FIGURE 1C — Volcano plot
#############################################################

# Ensure FDR exists (robustness)
if (!"adj.P.Val" %in% colnames(res2_ann)) {
  res2_ann$adj.P.Val <- p.adjust(res2_ann$p_value, method = "fdr")
}

# Ensure significance exists
if (!"significance" %in% colnames(res2_ann)) {
  res2_ann$significance <- dplyr::case_when(
    res2_ann$adj.P.Val < 0.05 & res2_ann$Estimate > 0  ~ "Upregulated",
    res2_ann$adj.P.Val < 0.05 & res2_ann$Estimate < 0 ~ "Downregulated",
    TRUE ~ "Not_significant"
  )
}

# -------------------------------
# Select top genes
# -------------------------------

DE_pos <- res2_ann %>%
  filter(adj.P.Val < 0.05, Estimate > 0) %>%
  arrange(desc(Estimate)) %>%
  slice_head(n = 20)

DE_neg <- res2_ann %>%
  filter(adj.P.Val < 0.05, Estimate < 0) %>%
  arrange(Estimate) %>%
  slice_head(n = 20)

genes <- bind_rows(DE_pos, DE_neg)$gene

# Labeling
res2_ann <- res2_ann %>%
  mutate(
    fig = ifelse(gene %in% genes, "TopGenes", "Other"),
    delabel = ifelse(fig == "TopGenes", gene, NA)
  )

# -------------------------------
# Threshold
# -------------------------------

fdr_threshold <- -log10(
  max(res2_ann$p_value[res2_ann$adj.P.Val <= 0.05], na.rm = TRUE)
)

# -------------------------------
# Plot
# -------------------------------

volplot <- ggplot(res2_ann, aes(x = Estimate, y = -log10(p_value))) +
  
  geom_point(
    aes(fill = significance),
    shape = 21,
    colour = "black",
    stroke = 0.3,
    size = 3,
    alpha = 0.9
  ) +
  
  geom_label_repel(
    data = dplyr::filter(res2_ann, !is.na(delabel)),
    aes(label = delabel),
    color = "black",
    max.overlaps = 50,
    force = 5,
    size = 4,
    show.legend = FALSE
  ) +
  
  scale_fill_manual(values = c(
    "Upregulated" = "#F77663",
    "Downregulated" = "#3D799E",
    "Not_significant" = "#f0f0f0"
  )) +
  
  geom_hline(yintercept = fdr_threshold,
             linetype = "dashed",
             color = "blue") +
  
  geom_vline(xintercept = 0,
             linetype = "dashed",
             color = "black") +
  
  theme_minimal() +
  theme(
    axis.title = element_text(size = 16),
    legend.position = "bottom",
    legend.text = element_text(size = 16),
    legend.title = element_blank()
  ) +
  
  labs(
    x = "Estimate",
    y = "-log10(p-value)"
  )

# -------------------------------
# Save
# -------------------------------

ggsave(
  filename = "figures/plots4figures/Figure1/Figure1C_volcano.png",
  plot = volplot,
  width = 8,
  height = 8,
  dpi = 300
)

#############################################################
# FIGURE 1D — Heatmap
#############################################################

# -------------------------------
# Meta
# -------------------------------

meta$Sample <- rownames(meta)

# Define new names
new_names <- c(
  "SSc_1","SSc_2","SSc_3","SSc_4","SSc_5","SSc_6","SSc_7","SSc_8",
  "CTRL_1","CTRL_2","CTRL_3","CTRL_4","CTRL_5","CTRL_6","CTRL_7","CTRL_8"
)

sample_names <- data.frame(
  og = unique(meta$sample),
  new = new_names
)

name_map <- setNames(sample_names$new, sample_names$og)

# -------------------------------
# Expression matrix
# -------------------------------

gene_plot <- genes

exp_heat <- exprs_norm[gene_plot, , drop = FALSE] %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame()

exp_heat$sample <- rownames(exp_heat)

# -------------------------------
# Average per sample
# -------------------------------

exp_ave <- aggregate(
  exp_heat[, 1:(ncol(exp_heat)-1)],
  list(exp_heat$sample),
  mean
)

colnames(exp_ave)[1] <- "sample"

# reshape
exp_ave <- as.data.frame(t(exp_ave))
colnames(exp_ave) <- exp_ave[1,]
exp_ave <- exp_ave[-1,]

# numeric matrix
exp_ave <- mutate_all(exp_ave, ~as.numeric(as.character(.)))
exp_ave <- as.matrix(exp_ave)

# -------------------------------
# Z-score
# -------------------------------

scale_rows <- function(x) t(scale(t(x)))

exp_ave_scale <- scale_rows(exp_ave)
exp_ave_scale[exp_ave_scale > 2] <- 2
exp_ave_scale[exp_ave_scale < -2] <- -2

# -------------------------------
# Rename samples
# -------------------------------

colnames(exp_ave_scale) <- ifelse(
  colnames(exp_ave_scale) %in% names(name_map),
  name_map[colnames(exp_ave_scale)],
  colnames(exp_ave_scale)
)

# -------------------------------
# Order samples: CTRL → SSc
# -------------------------------

ordered_cols <- colnames(exp_ave_scale)[
  order(
    grepl("^SSc", colnames(exp_ave_scale)),   # CTRL first
    as.numeric(gsub(".*_", "", colnames(exp_ave_scale)))
  )
]

exp_ave_scale <- exp_ave_scale[, ordered_cols]

# -------------------------------
# Plot and Save
# -------------------------------

pheatmap(
  exp_ave_scale,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_col = 10,
  angle_col = 45,
  treeheight_row = 0,
  treeheight_col = 0,
  fontsize = 12,
  border_color = NA,
  filename = "figures/plots4figures/Figure1/Figure1D.png",
  width = 9,
  height = 10
)


#############################################################
# Figure 1E - GSEA (Pseudobulk)
#############################################################


hallmark_pseub <- read.csv("ShinyCD4/data/GSEA_TFactivity/2026-01-31_SSc_CD4_GSEA_pseudobulk_Hallmark.csv")
reactome_pseub <- read.csv("ShinyCD4/data/GSEA_TFactivity/2026-01-31_SSc_CD4_GSEA_pseudobulk_C2.csv")

# Merge datasets
gsea_res <- bind_rows(hallmark_pseub, reactome_pseub)

# -------------------------------
# Basic sanity checks
# -------------------------------

stopifnot(all(c("pathway", "NES", "padj") %in% colnames(gsea_res)))

# -------------------------------
# Filter significant pathways
# -------------------------------

gsea_filtered <- gsea_res %>%
  filter(padj < 0.1)

# -------------------------------
# Select top enriched pathways
# -------------------------------

# Top positively enriched pathways (SSc upregulated)
top_pos <- gsea_filtered %>%
  filter(NES > 0) %>%
  arrange(desc(NES)) %>%
  slice_head(n = 20)

# Top negatively enriched pathways (CTRL upregulated)
top_neg <- gsea_filtered %>%
  filter(NES < 0) %>%
  arrange(NES) %>%
  slice_head(n = 20)

# Combine results
top_combined <- bind_rows(top_pos, top_neg)

# -------------------------------
# Improve pathway readability
# -------------------------------

top_combined <- top_combined %>%
  mutate(
    pathway_clean = pathway %>%
      str_remove("^HALLMARK_") %>%
      str_remove("^REACTOME_") %>%
      str_replace_all("_", " ") %>%
      str_to_title()
  )

# -------------------------------
# Generate GSEA barplot
# -------------------------------

gsea_plot <- ggplot(top_combined, aes(x = reorder(pathway_clean, NES), y = NES)) +
  geom_col(
    aes(fill = NES > 0),
    alpha = 0.85,
    width = 0.75,
    color = "black"
  ) +
  scale_fill_manual(values = c("NES > 0" = "#FF7F50", "NES < 0" = "#6495ED"), 
                    name = "FDR 10%") +
  coord_flip() +
  labs(
    x = "",
    y = "Normalized Enrichment Score (NES)",
    title = "GSEA top pathways in SSc vs CTRL"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 10, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.position = "top",
    legend.text = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )


# Prepare data for plotting
gsea_dotplot_df <- top_combined %>%
  mutate(
    pathway_clean = pathway %>%
      str_remove("^HALLMARK_") %>%
      str_remove("^REACTOME_") %>%
      str_replace_all("_", " ") %>%
      str_to_title(),
    log10_padj = -log10(padj)
  )

# Dotplot
gsea_dotplot <- ggplot(
  gsea_dotplot_df,
  aes(
    x = NES,
    y = reorder(pathway_clean, NES)
  )
) +
  geom_point(
    aes(
      size = log10_padj,
      color = NES
    )
  ) +
  scale_color_gradient2(
    low = "#6495ED",   # blue (CTRL enriched)
    mid = "white",
    high = "#FF7F50",  # red (SSc enriched)
    midpoint = 0,
    name = "NES"
  ) +
  scale_size_continuous(
    name = "-log10(FDR)"
  ) +
  labs(
    x = "Normalized Enrichment Score (NES)",
    y = "",
    title = "GSEA top pathways in SSc vs CTRL"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 10, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.position = "right",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# -------------------------------
# Save figure
# -------------------------------

ggsave(
  filename = "figures/plots4figures/Figure1/Figure1E.png",
  plot = gsea_plot,
  width = 12,
  height = 10,
  dpi = 300
)


#############################################################
# Figure 1F - TF activity (Pseudobulk)
#############################################################

f_contrast_acts <- read.csv("ShinyCD4/data/GSEA_TFactivity/2026-01-31_SSc_CD4_GSEA_pseudobulk_TFact.csv")

colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")[c(2, 10)])

p <- ggplot(data = f_contrast_acts, 
            mapping = aes(x = stats::reorder(source, score), 
                          y = score)) + 
  geom_bar(mapping = aes(fill = score),
           color = "black",
           stat = "identity") +
  scale_fill_gradient2(low = colors[1], 
                       mid = "whitesmoke", 
                       high = colors[2], 
                       midpoint = 0) + 
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 16),
        axis.text.x = element_text(angle = 45, 
                                   hjust = 1, 
                                   size = 14, 
                                   face = "bold"),
        axis.text.y = element_text(size = 13, face = "bold"),
        legend.title = element_text(size = 16,face = "bold"),
        legend.position = "top",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  lab(title = "Transcription Factor activity score",
    x= "TFs",
    y= "score")

p

ggsave(p,filename = "figures/plots4figures/Figure1/Figure1F.png",
       width = 8,
       height = 6, 
       dpi = 300)

#############################################################
# Figure 1G - Modules score
#############################################################

# -------------------------------
# Helper function: module volcano
# -------------------------------

moduleVolcano <- function(res_df, genes, label_name, color_val, fdr_threshold) {
  
  res_df$group <- ifelse(res_df$gene %in% genes, label_name, "Other")
  res_df$delabel <- ifelse(res_df$group != "Other", res_df$gene, NA)
  
  p <- ggplot(res_df, aes(x = Estimate, y = -log10(p_value))) +
    
    geom_point(
      data = subset(res_df, group == "Other"),
      shape = 21,
      fill = "#f0f0f0",
      colour = "black",
      stroke = 0.2,
      size = 2,
      alpha = 0.6
    ) +
    
    geom_point(
      data = subset(res_df, group != "Other"),
      aes(fill = group),
      shape = 21,
      colour = "black",
      stroke = 0.4,
      size = 3.5
    ) +
    
    geom_label_repel(
      data = subset(res_df, gene %in% genes),
      aes(label = gene),
      size = 3.5,
      max.overlaps = Inf,
      force = 5,
      show.legend = FALSE
    ) +
    
    scale_fill_manual(values = setNames(color_val, label_name)) +
    
    geom_hline(yintercept = fdr_threshold, linetype = "dashed", color = "blue") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    
    theme_minimal() +
    theme(legend.position = "") +
    
    labs(
      x = "Estimate",
      y = "-log10(p-value)"
    )
  
  return(p)
}

gene_score <- function(genes, meta, exprs_norm, score_name){
  
  ## -------------------------------
  ## 1. Filtering the genes
  ## -------------------------------
  genes_present <- intersect(genes, rownames(exprs_norm))
  genes_missing <- setdiff(genes, genes_present)
  
  message("Genes present: ", length(genes_present))
  if (length(genes_missing) > 0) {
    message("Genes missing: ", paste(genes_missing, collapse = ", "))
  }
  
  if (length(genes_present) < 2) {
    stop("Not enough genes present to compute a score")
  }
  
  ## -------------------------------
  ## 2. Subset expression
  ## -------------------------------
  exp_heat <- exprs_norm[genes_present, , drop = FALSE]
  exp_heat <- as.matrix(exp_heat)
  exp_heat <- t(exp_heat)
  exp_heat <- as.data.frame(exp_heat)
  
  exp_heat$sample <- as.character(meta$Sample)
  
  ## -------------------------------
  ## 3. Average per sample
  ## -------------------------------
  exp_ave <- aggregate(
    exp_heat[, 1:(ncol(exp_heat)-1)], 
    list(exp_heat$sample), 
    mean
  )
  
  colnames(exp_ave)[1] <- "sample"
  exp_ave <- as.data.frame(t(exp_ave))
  
  colnames(exp_ave) <- as.character(t(exp_ave[1,]))
  exp_ave <- exp_ave[-1,]
  row_names <- rownames(exp_ave)
  
  exp_ave <- mutate_all(exp_ave, function(x) as.numeric(as.character(x)))
  
  rownames(exp_ave) <- row_names
  exp_ave <- as.matrix(exp_ave)
  
  ## -------------------------------
  ## 5. Row-wise Z-score
  ## -------------------------------                      
  
  scale_rows <- function(x) t(scale(t(x)))
  exp_ave_scale <- scale_rows(exp_ave) # Z-score
  exp_ave_scale[exp_ave_scale > 2] <- 2
  exp_ave_scale[exp_ave_scale < -2] <- -2
  
  
  score <- colMeans(exp_ave_scale[genes, , drop = FALSE])                          
  
  return(score)
}

moduleBoxPlot <- function(meta,score,title){
  pltitle <- paste0(title," score")
  module_Plot <- ggplot(meta, aes(x=factor(condition, levels=c("CTRL","SSc")), y=score, fill=condition)) + 
    geom_boxplot(width=0.5, alpha = 0.75)+
    geom_point(alpha = 0.9, 
               shape=21, 
               colour="black", 
               stroke=0.3, 
               size=2.5) +
    labs(title = pltitle, 
         x = "", 
         y = "z-score",
         fill="") +
    theme_bw()+
    theme(legend.position = '',
          axis.text.x = element_text(size=12,face = "bold"),
          axis.text.y = element_text(size=12),
          axis.title.y = element_text(size=20),
          title = element_text(size = 18),
          plot.margin = unit(c(1, 1, 1, 1), "mm"))+
    scale_fill_manual(guide = guide_legend(reverse=TRUE), 
                      values=c('#D8BFD8','#C3E5C3'))+
    stat_compare_means(comparisons = list(c("CTRL","SSc")) , label = "p.signif", method = "wilcox.test")
  return(module_Plot)
}

# -------------------------------
# FDR threshold
# -------------------------------

fdr_threshold <- -log10(max(res2_ann$p_value[res2_ann$adj.P.Val <= 0.05], na.rm = TRUE))

# -------------------------------
# 1. Exhaustion module
# -------------------------------

genes_exh <- c("PDCD1","HAVCR2","LAG3","TOX","TOX2","NR4A1","NR4A2","NR4A3","PRDM1","ENTPD1")

score_exh <- gene_score(genes_exh, meta, exprs_norm, "Exhaustion")

df_exh <- data.frame(Sample = names(score_exh), score = score_exh) %>%
  left_join(meta_ssc, by = "Sample")

plot_exh_box <- moduleBoxPlot(df_exh, df_exh$score, "Exhaustion")

plot_exh_vol <- moduleVolcano(res2_ann, genes_exh, "Exhaustion", "#31a354", fdr_threshold)

# -------------------------------
# 2. Activation module
# -------------------------------

genes_act <- c("HLA-DRA","HLA-DRB1","HLA-DQA1","HLA-DQB1","HLA-DPA1","HLA-DPB1","CD69","IL2RA")

score_act <- gene_score(genes_act, meta_ssc, exprs_norm, "Activation")

df_act <- data.frame(Sample = names(score_act), score = score_act) %>%
  left_join(meta_ssc, by = "Sample")

plot_act_box <- moduleBoxPlot(df_act, df_act$score, "Activation")

plot_act_vol <- moduleVolcano(res2_ann, genes_act, "Activation", "#3182bd", fdr_threshold)

# -------------------------------
# 3. IFN module
# -------------------------------

IFNa_genes <- gsea.H %>%
  filter(pathway == "HALLMARK_INTERFERON_ALPHA_RESPONSE") %>%
  pull(leadingEdge) %>%
  unlist()

IFNg_genes <- gsea.H %>%
  filter(pathway == "HALLMARK_INTERFERON_GAMMA_RESPONSE") %>%
  pull(leadingEdge) %>%
  unlist()

genes_ifn <- unique(c(IFNa_genes, IFNg_genes))

score_ifn <- gene_score(genes_ifn, meta_ssc, exprs_norm, "IFN")

df_ifn <- data.frame(Sample = names(score_ifn), score = score_ifn) %>%
  left_join(meta_ssc, by = "Sample")

plot_ifn_box <- moduleBoxPlot(df_ifn, df_ifn$score, "IFN")

plot_ifn_vol <- moduleVolcano(res2_ann, genes_ifn, "IFN", "#c51b8a", fdr_threshold)

# -------------------------------
# Combine plots
# -------------------------------

box_row <- plot_exh_box + plot_act_box + plot_ifn_box
vol_row <- plot_exh_vol + plot_act_vol + plot_ifn_vol

figure1G <- box_row / vol_row

# -------------------------------
# Save figure
# -------------------------------

ggsave(
  filename = "figures/plots4figures/Figure1/Figure1G.png",
  plot = figure1G,
  width = 14,
  height = 8,
  dpi = 300
)
#!/usr/bin/env Rscript

# ------------------------------------------------------------
# Figure 2 generation script
# CD4 T cell single-cell RNA-seq analysis in SSc
# ------------------------------------------------------------

# Load required libraries quietly (suppress startup messages)
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
  library(readr)
  library(stringr)
  library(patchwork)
  library(fgsea)
  library(msigdbr)
  library(cowplot)
  library(ggh4x)
  library(viridis)
  library(ggpubr)
})

# ------------------------------------------------------------
# Create output directory
# ------------------------------------------------------------

# Create directory to store Figure 2 outputs (recursive to ensure full path exists)
dir.create("figures/plots4figures/Figure2", recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Load metadata
# ------------------------------------------------------------

# Load main metadata (cells x annotations)
meta <- read.csv("data/20260124_SSc_CD4_meta_v3.csv",
                 row.names = 1,
                 stringsAsFactors = FALSE)

# Define cluster labels from R09 annotation
meta$clusters <- as.character(meta$R09)

# Load normalized expression matrix (genes x cells)
exprs_norm <- readRDS("ShinyCD4/data/20260124_SSc_CD4_counts_norm.rds")

# ------------------------------------------------------------
# Color palette
# ------------------------------------------------------------

# Define color palette for clusters and conditions
custom_colors <- c(
  "0" = "#CDB4DB", 
  "1" = "#F78181", 
  "2" = "#82CFFD",
  "3" = "#8FCFD1", 
  "4" = "#AEB6BF", 
  "5" = "#F7ACB4", 
  "6" = "#FFB347",
  "7" = "#7FFFD4", 
  "8" = "#F4D03F", 
  "9" = "#F8D5E1",
  "10" = "#FFA07A",
  "11" = "#2ECC71", 
  "CTRL" = "#D8BFD8", 
  "SSc" = "#C3E5C3"
)

# Separate cluster color palette
clusterCols <- c(
  "0" = "#CDB4DB", "1" = "#F78181", "2" = "#82CFFD","3" = "#8FCFD1", "4" = "#AEB6BF", 
  "5" = "#F7ACB4", "6" = "#FFB347","7" = "#7FFFD4", "8" = "#F4D03F", "9" = "#F8D5E1",
  "10" = "#FFA07A","11" = "#2ECC71")

# Cluster labels including representative marker genes
clusterNum <- c('0 CCR7/EEF1B2/TUBA1A', '1 MALAT1/EEF1B2/CCR7', '2 FOSB/FOS/IER2','3 S100A4/KLRB1/LTB','4 YPEL5/EIF1/SARAF', '5 VIM/ANXA1/S100A4','6 FTH1/JUNB/SRGN','7 SRGN/GZMK/ZFP36','8 TPT1/DDX5/EEF2','9 S100A4/LMNA/S100A10', '10 FOXP3/RTKN2/TIGIT',  '11 NKG7/CCL5/GZMH')

# ------------------------------------------------------------
# Layer 1 - Clustering and abundance
# ------------------------------------------------------------

# Abundance ----

# Define cluster ordering
cluster_order <- c(0:11)

# Compute cluster proportions per sample
abundance_sample <- meta %>%
  count(Sample, condition, R09) %>%
  group_by(Sample) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# Compute cluster proportions per condition
abundance_condition <- meta %>%
  dplyr::count(condition, R09) %>%
  group_by(condition) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# Assign human-readable cluster names (cluster + marker genes)
abundance_sample <- abundance_sample %>% 
  mutate(cnames = case_when(
    R09==0 ~ '0 CCR7/EEF1B2/TUBA1A', 
    R09==1 ~ '1 MALAT1/EEF1B2/CCR7', 
    R09==2 ~ '2 FOSB/FOS/IER2',
    R09==3 ~ '3 S100A4/KLRB1/LTB',
    R09==4 ~ '4 YPEL5/EIF1/SARAF', 
    R09==5 ~ '5 VIM/ANXA1/S100A4',
    R09==6 ~ '6 FTH1/JUNB/SRGN',
    R09==7 ~ '7 SRGN/GZMK/ZFP36',
    R09==8 ~ '8 TPT1/DDX5/EEF2',
    R09==9 ~ '9 S100A4/LMNA/S100A10',
    R09==10 ~ '10 FOXP3/RTKN2/TIGIT',
    R09==11 ~ '11 NKG7/CCL5/GZMH'
  ))

# ------------------------------------------------------------
# UMAP annotation + Bar plot(Figure 2A)
# ------------------------------------------------------------

# Compute cluster centroids in UMAP space
centroids <- meta %>%
  group_by(clusters) %>%
  summarize(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2))

# UMAP visualization of clusters
umap_clusters <- ggplot(meta, aes(x = UMAP_1, y = UMAP_2, color = clusters)) +
  geom_point(color = "black", size = 0.9) +
  geom_point(aes(color = clusters), size = 0.7) +
  scale_color_manual(values = clusterCols) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_blank(),
    legend.text = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = "black", size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 14),
    aspect.ratio = 1
  ) +
  ggtitle("") +
  geom_text(data = centroids, 
            aes(x = UMAP_1, y = UMAP_2, label = clusters), 
            color = "black", 
            size = 5, 
            fontface = "bold")

# Barplot of cluster composition per sample
ab_sample_clusters <- ggplot(abundance_sample, 
                             aes(fill=factor(R09, levels =rev(cluster_order)), x=Sample, y=prop)
) +
  geom_bar(stat = "identity", color='black')+
  labs(fill="",x=NULL,y="Fraction of cells")+
  theme_void()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.ticks = element_blank(),
        axis.text=element_text(size=11,
                               face="bold"),
        axis.text.x = element_text(angle = 45, 
                                   vjust = 0.9, 
                                   hjust=1),
        axis.title=element_blank())+
  scale_fill_manual(guide = guide_legend(reverse=TRUE), 
                    values=clusterCols)

# Legend ----

# Create data for custom legend (cluster colors + labels)
data <- data.frame(
  y = 11:0,
  color = clusterCols,
  number = clusterNum
)

# Custom legend plot
legendClusters <- ggplot(data, aes(x = 1, y = y, fill = color)) +
  geom_tile(width = 0.6, height = 0.6, color = "black") +
  geom_text(aes(x = 1.5, label = number), hjust = 0, color = "black", size = 4, fontface = "bold") +
  scale_fill_identity() +
  coord_fixed(ratio = 1) +
  scale_x_continuous(expand = c(0, 0), breaks = NULL, limits = c(0.5, 12)) +
  scale_y_continuous(expand = c(0, 0), breaks = NULL) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.margin = unit(c(1, 1, 1, 1), "mm")
  )

# Combine UMAP, barplot and legend
comb.1 <- wrap_plots(umap_clusters, ab_sample_clusters, legendClusters)

# Save Figure 2A
ggsave(comb.1, 
       filename = "figures/plots4figures/Figure2/figure2A.png",
       width = 7,
       height = 5,
       dpi = 300)

# ------------------------------------------------------------
# Boxplot (Figure 2B)
# ------------------------------------------------------------

# Boxplot of selected clusters (3, 4, 9) comparing conditions
fig2B <- ggplot(abundance_sample %>% filter(R09 %in% c("3","4","9")), 
                aes(x = factor(condition, levels=c("CTRL","SSc")), 
                    y = prop, 
                    fill = condition)) +
  geom_boxplot(width = 0.5, alpha = 0.75) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.6) +
  stat_compare_means(comparisons = list(c("CTRL","SSc")) , label = "p.signif", method = "wilcox.test")+
  facet_wrap2(
    ~ factor(cnames, levels = c("3 S100A4/KLRB1/LTB","4 YPEL5/EIF1/SARAF","9 S100A4/LMNA/S100A10")),
    scales = "free",
    ncol = 3,
    strip = strip_themed(
      background_x = elem_list_rect(
        fill = c("#FFB347", "#AEB6BF","#F8D5E1"),
        colour = "black"
      )
    )
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(size = 10),
    axis.title = element_blank(),
    plot.margin = unit(c(1, 1, 1, 1), "mm"),
    strip.text = element_text(face = "bold", size = 12)
  ) +
  scale_fill_manual(
    values = c("#D8BFD8", "#C3E5C3")
  ) + 
  ylim(0,0.20)

# Save Figure 2B
ggsave(plot = fig2B, 
       filename = "figures/plots4figures/Figure2/figure2B.png", 
       width = 7.4,
       height = 2,
       dpi = 300)

# ------------------------------------------------------------
# DEA quantification (Figure 2D)
# ------------------------------------------------------------

# Helper functions to clean DEA filenames
clean_dea_name <- function(filename, prefix) {
  basename <- tools::file_path_sans_ext(basename(filename))
  cleaned <- sub("scRNAseq_SScCD4_", "", basename)
  cleaned <- sub("_DEG_R09_byCON_1500HVG_20PC_it92", "", cleaned)
}

clean_dea_CT_name <- function(filename, prefix) {
  basename <- tools::file_path_sans_ext(basename(filename))
  cleaned <- sub("scRNAseq_SScCD4_", "", basename)
  cleaned <- sub("_Soskic_DEG_R09_byCON_1500HVG_20PC_it92", "", cleaned)
}

# Convert Scanpy output to DEA table
scanpy_to_dea <- function(df,cluster_id) {
  out <- data.frame(
    gene   = df$SSc_names,
    logFC  = df$SSc_logfoldchanges,
    pval   = df$SSc_pvals,
    padj   = df$SSc_pvals_adj,
    zscore   = df$SSc_scores,
    cluster = cluster_id,
    stringsAsFactors = FALSE
  )
  
  out <- out[!is.na(out$gene), ]
  return(out)
}

# Load DEA files
base_path <- "ShinyCD4/data/dea_sc/"

cluster_files <- list.files(file.path(base_path, "cluster"), pattern = "*.csv", full.names = TRUE)
cluster_list <- lapply(cluster_files, read.csv)
names(cluster_list) <- vapply(cluster_files, clean_dea_name, character(1))

# Build DEA table across clusters
cluster_id <- names(cluster_list)

dea_clusters <- lapply(names(cluster_list), function(cl) {
  scanpy_to_dea(cluster_list[[cl]],cl)
})

names(dea_clusters) <- names(cluster_list)

dea_cluster_df <- do.call(rbind,dea_clusters)

cluster_order <- paste0("Cluster",c(0:11))

# Filter significant genes and classify direction
dea_sig <- dea_cluster_df %>%
  filter(padj < 0.05) %>%
  mutate(direction = ifelse(logFC > 0, "Up", "Down"))

# Count DE genes per cluster and direction
dea_counts <- dea_sig %>%
  count(cluster, direction) %>%
  mutate(
    n_signed = ifelse(direction == "Up", n, -n)
  )

dea_counts$cluster <- factor(dea_counts$cluster, levels = rev(cluster_order))

# Plot DEA counts
fig2D <- ggplot(dea_counts, aes(x = cluster, y = n_signed, fill = direction)) +
  geom_col(width = 0.7, color = "black") +
  scale_fill_manual(
    values = c("Up" = "#F77663", "Down" = "#3D799E"),
    breaks = c("Up", "Down")
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 11),
    axis.text.x = element_text(angle = 45,vjust = 0.6),
    plot.title = element_text(hjust = .3,size = 14),
    legend.position = "top"
  ) +
  labs(
    y = "Number of genes",
    x = "",
    fill = "Direction",
    title = "Significant DE genes \n per cluster in SSc"
  ) + 
  coord_flip()

# Save Figure 2D
ggsave(plot = fig2D, 
       filename = "figures/plots4figures/Figure2/figure2D.png",
       width = 2.8,
       height = 4.4,
       dpi = 300)

# ------------------------------------------------------------
# Hallmark pathways (Figure 2E)
# ------------------------------------------------------------

# Load GSEA results and assign names
gsea_res <- list.files("results/GSEA_scanRes/", pattern = "\\.csv$", full.names = TRUE) %>%
  setNames(gsub("2026-02-03_SSc_CD4_GSEA_|\\.csv$", "", basename(.))) %>%
  lapply(read_csv)

g_files <- list.files(file.path("results/GSEA_scanRes/"), pattern = "*.csv", full.names = TRUE)

# Clean GSEA file names
clean_gsea_name <- function(filename, prefix) {
  basename <- tools::file_path_sans_ext(basename(filename))
  cleaned <- sub("^2026-02-03_SSc_CD4_GSEA_", "", basename)
}

names(gsea_res) <- vapply(g_files, clean_gsea_name, character(1))

gsea_res_df <- do.call(rbind, gsea_res)

cluster_names <- paste0("cluster_",c(0:11))

# Load pseudobulk hallmark results
hallmark_pseub <- read.csv("ShinyCD4/data/GSEA_TFactivity/2026-01-31_SSc_CD4_GSEA_pseudobulk_Hallmark.csv")

# Select top positive and negative pathways
hallmark_pseub_pos <- hallmark_pseub %>%
  filter(padj < 0.01, NES > 0) %>%
  arrange(NES) %>%
  slice_head(n=6)

hallmark_pseub_neg <- hallmark_pseub %>%
  filter(padj < 0.01, NES < 0) %>%
  arrange(NES) %>%
  slice_head(n=6)

hallmark_pseub_fil <- rbind(hallmark_pseub_pos,hallmark_pseub_neg)

hallmark_pseub_v <- as.vector(hallmark_pseub_fil$pathway)

# Filter cluster-level hallmark pathways
gsea_hallmark <- gsea_res_df %>%
  filter(pathway %in% hallmark_pseub_v) %>%
  select(-...1) %>% 
  filter(grepl("^HALLMARK", pathway)) %>%
  filter(padj < 0.1)

# Convert to wide format
gsea_hallmark_wide <- gsea_hallmark %>%
  filter(padj < 0.1)%>%
  select(pathway,NES,cluster) %>%
  pivot_wider(
    names_from = cluster,
    values_from = NES
  )

cluster4path <- c("Cluster0","Cluster1","Cluster2","Cluster3","Cluster4","Cluster5",
                  "Cluster6","Cluster7","Cluster8","Cluster9","Cluster10","Cluster11")

gsea_hallmark_wide <- gsea_hallmark_wide[, c("pathway",0:11)]

colnames(gsea_hallmark_wide) <- c("pathway",cluster4path)

# Select top pathways with least missing values
top_pathways_h <- gsea_hallmark_wide %>%
  ungroup() %>%
  mutate(n_na = rowSums(is.na(across(where(is.numeric))))) %>%
  arrange(n_na) %>%
  slice_head(n = 20)

# Clean pathway names
top_pathways_h$pathway <- sub("^HALLMARK_", "", top_pathways_h$pathway)

# Convert to long format for plotting
top_path_longer_h <- top_pathways_h %>%
  pivot_longer(
    cols = cluster4path,
    names_to = "cluster",
    values_to = "NES"
  ) %>%
  mutate(
    pathway = factor(
      pathway,
      levels = top_pathways_h %>%
        arrange(desc(n_na)) %>%
        pull(pathway)
    )
  )

# Dotplot of hallmark pathways
fig2E <- ggplot(
  top_path_longer_h %>% filter(!is.na(NES)),
  aes(x = factor(cluster, levels = cluster4path), y = pathway)
) +
  geom_point(
    aes(
      size = abs(NES),
      fill = ifelse(NES > 0, "NES > 0", "NES < 0")
    ),
    shape = 21,
    stroke = 0.5
  ) +
  theme_bw() + 
  theme(
    axis.title.y = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = .6,size = 14),
    axis.text.x = element_text(size = 9, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 9, face = "bold")
  ) +
  labs(
    x = "",
    y = "",
    title = "Hallmark Pathways across Clusters",
    size = "abs(NES)",
    fill = "Direction"
  ) +
  scale_fill_manual(
    values = c("NES > 0" = "#F77663", "NES < 0" = "#3D799E")
  ) +
  scale_size_continuous(name = "NES - FDR 10%")

# Save Figure 2E
ggsave(plot=fig2E, 
       filename = "figures/plots4figures/Figure2/figure2E.png",
       width = 7,
       height = 2.9,
       dpi = 300)

# ------------------------------------------------------------
# IFN alpha and gamma heatmap (Figure 2F)
# ------------------------------------------------------------

# Define IFN-alpha and IFN-gamma gene signatures
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

# Keep genes present in dataset
IFNa_present <- IFNa_genes[IFNa_genes %in% rownames(exprs_norm)]
IFNg_present <- IFNg_genes[IFNg_genes %in% rownames(exprs_norm)]

# Compute per-cell signature scores
meta$IFNa_score <- Matrix::colMeans(exprs_norm[IFNa_present, , drop = FALSE])
meta$IFNg_score <- Matrix::colMeans(exprs_norm[IFNg_present, , drop = FALSE])

# Aggregate scores per cluster and condition
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

# Order clusters
score_wide$clusters <- factor(score_wide$clusters,
                              levels = as.character(0:11))

score_wide <- score_wide %>%
  arrange(clusters)

# Build matrix for plotting
score_mat <- as.matrix(score_wide[, -1])
rownames(score_mat) <- paste0("cluster_", score_wide$clusters)

# Order columns explicitly
score_mat <- score_mat[, c(
  "IFNa_CTRL","IFNa_SSc",
  "IFNg_CTRL","IFNg_SSc"
)]

# Annotation for heatmap
annotation_col <- data.frame(
  Signature = c("IFNa","IFNa","IFNg","IFNg"),
  Condition = c("CTRL","SSc","CTRL","SSc")
)

rownames(annotation_col) <- colnames(score_mat)

# Annotation colors
annotation_colors <- list(
  Signature = c(IFNa = "#5E548E", IFNg = "#1B998B"),
  Condition = c(CTRL = "#D8BFD8", SSc = "#C3E5C3")
)

score_mat_plot <- score_mat

# Ensure consistent column naming
colnames(score_mat_plot) <- c(
  "IFNa_CTRL","IFNa_SSc",
  "IFNg_CTRL","IFNg_SSc"
)

# Transpose matrix for plotting
score_mat_t <- t(score_mat_plot)

# Align annotation with rows
annotation_row <- annotation_col
rownames(annotation_row) <- rownames(score_mat_t)

# Define gap between IFNa and IFNg
gaps_row <- 2  

# Heatmap (without legend)
ph <- pheatmap(
  t(score_mat_plot),
  color = viridis(100),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  border_color = "black",
  show_rownames = FALSE,
  
  annotation_row = annotation_col,
  annotation_colors = annotation_colors,
  
  gaps_row = gaps_row,
  
  angle_col = 45,
  fontsize_col = 9,
  
  legend = FALSE,
  
  annotation_legend = TRUE,
  
  legend_labels = NULL,
  
  fontsize = 10,
  
  annotation_names_row = FALSE,
  
  annotation_names_col = FALSE,
  main = "Cluster-level IFN-α and IFN-γ signature score"
)

# Extract range for legend
min_val <- min(score_mat_plot)
max_val <- max(score_mat_plot)

# Create legend gradient
legend_df <- data.frame(
  x = 1,
  y = seq(min_val, max_val, length.out = 100)
)

legend_plot <- ggplot(legend_df, aes(x = x, y = y, fill = y)) +
  geom_tile() +
  scale_fill_viridis_c(
    name = "Mean of expression in group",
    limits = c(min_val, max_val),
    guide = guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0.5,
      barwidth = unit(5, "cm"),
      barheight = unit(0.4, "cm")
    )
  ) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )

# Extract legend only
legend_only <- cowplot::get_legend(legend_plot)

# Combine heatmap and legend
final_plot <- cowplot::plot_grid(
  ph$gtable,
  legend_only,
  ncol = 1,
  rel_heights = c(4, 1)
)

# Display final plot
final_plot

# Save Figure 2F
ggsave(plot = final_plot, 
       filename = "figures/plots4figures/Figure2/figure2F.png", 
       width = 7.2,
       height = 3,
       dpi = 300)
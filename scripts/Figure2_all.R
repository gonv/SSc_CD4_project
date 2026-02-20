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
  library(readr)
  library(stringr)
  library(patchwork)
  library(fgsea)
  library(msigdbr)
  library(cowplot)
  library(ggh4x)
})

# Loading files----

meta <- read.csv("data/20260124_SSc_CD4_meta_v3.csv", row.names = 1)
meta$clusters <- as.character(meta$R09)

# Colors and vectors----
cluster_order <- c(0:11)

custom_colors <- c(
  "0" = "#CDB4DB", "1" = "#F78181", "2" = "#82CFFD", "6" = "#FFB347", "4" = "#AEB6BF",
  "5" = "#F7ACB4", "7" = "#7FFFD4", "10" = "#FFA07A", "8" = "#F4D03F", "9" = "#F8D5E1",
  "3" = "#8FCFD1", "11" = "#2ECC71", "CTRL" = "#D8BFD8", "SSc" = "#C3E5C3",
  "Tn"= "#a9dfbb", "Tcm"= "#6cd0c4", "Treg"= "#fc6daf", "Tem"= "#3faccd",
  "Temra"= "#ff0000", "Th17"= '#d7ff55', "Th2"= '#FF9F55', "Thf"= '#559FFF',
  "Th1_GZMK"= '#55FF9F', "Th1_GZMB"= '#FFC255',"Th22"='#55FF9F',"Th1"='#FFC255'
)

clusterCols <- c(
  "0" = "#CDB4DB", "1" = "#F78181", "2" = "#82CFFD", "6" = "#FFB347", "4" = "#AEB6BF",
  "5" = "#F7ACB4", "7" = "#7FFFD4", "10" = "#FFA07A", "8" = "#F4D03F", "9" = "#F8D5E1",
  "3" = "#8FCFD1", "11" = "#2ECC71")

clusterNum <- c('0_CCR7_EEF1B2_TUBA1A', '1_MALAT1_EEF1B2_CCR7', '2_FOSB_FOS_IER2','3_S100A4_KLRB1_LTB','4_YPEL5_EIF1_SARAF', '5_VIM_ANXA1_S100A4','6_FTH1_JUNB_SRGN','7_SRGN_GZMK_ZFP36','8_TPT1_DDX5_EEF2','9_S100A4_LMNA_S100A10', '10_FOXP3_RTKN2_TIGIT',  '11_NKG7_CCL5_GZMH')

# Abundance ----
abundance_sample <- meta %>%
  count(Sample, condition, R09) %>%
  group_by(Sample) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

abundance_condition <- meta %>%
  dplyr::count(condition, R09) %>%
  group_by(condition) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# Plots ----
## Box plots ----
ggplot(abundance_sample, aes(x=factor(condition, levels=c("CTRL","SSc")), y=prop, fill=condition)) + 
  facet_wrap(~R09, ncol = 5)+
  geom_boxplot(width=0.5, alpha = 0.75)+
  geom_dotplot(binaxis='y', stackdir='center',binwidth = 0.0035)+
  theme_bw()+
  theme(legend.position = '',
        axis.text.x = element_text(size=10,face = "bold"),
        axis.text.y = element_text(size=10),
        axis.title = element_blank(),
        title = element_text(size = 10),
        plot.margin = unit(c(1, 1, 1, 1), "mm"))+
  scale_fill_manual(guide = guide_legend(reverse=TRUE), 
                    values=custom_colors)

ggplot(abundance_sample, aes(x = factor(condition, levels=c("CTRL","SSc")), y = prop, fill = condition)) +
  geom_boxplot(width = 0.5, alpha = 0.75) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.6) +
  facet_wrap2(
    ~ factor(R09, levels = cluster_order),
    scales = "free",
    nrow = 4,
    strip = strip_themed(
      background_x = elem_list_rect(
        fill = c("#CDB4DB","#F78181", "#82CFFD", "#FFB347", "#AEB6BF",
                 "#F7ACB4", "#7FFFD4", "#FFA07A", "#F4D03F", "#F8D5E1",
                 "#8FCFD1", "#2ECC71"),
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
    plot.margin = unit(c(1, 1, 1, 1), "mm")
  ) +
  scale_fill_manual(
    values = c("#D8BFD8", "#C3E5C3")
  )
#
## Bar plots ----
# By condition
ggplot(abundance_condition, aes(fill=factor(R09, levels =cluster_order), x=condition, y=prop)) +
  geom_bar(stat = "identity", color='black')+
  labs(fill="Cell_Type",x=NULL,y="Fraction of cells")+
  theme_classic()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #legend.key = element_rect(fill = "lightblue", color = NA),
        legend.key.size = unit(0.5, "cm"),
        legend.key.width = unit(0.5,"cm"), 
        axis.ticks = element_blank(),
        axis.text=element_text(size=14,face="bold"),
        axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1),
        axis.title=element_blank())+
  scale_fill_manual(guide = guide_legend(reverse=TRUE), 
                    values=custom_colors)
# by sample
ab_sample_clusters <- ggplot(abundance_sample, aes(fill=factor(R09, levels =rev(cluster_order)), x=Sample, y=prop)) +
  geom_bar(stat = "identity", color='black')+
  labs(fill="Cell_Type",x=NULL,y="Fraction of cells")+
  theme_void()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.ticks = element_blank(),
        axis.text=element_text(size=14,face="bold"),
        axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1),
        axis.title=element_blank())+
  scale_fill_manual(guide = guide_legend(reverse=TRUE), 
                    values=custom_colors)


# UMAPs ----
centroids <- meta %>%
  group_by(clusters) %>%
  summarize(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2))

umap_clusters <- ggplot(meta, aes(x = UMAP_1, y = UMAP_2, color = clusters)) +
  geom_point(color = "black", size = 0.9) +
  geom_point(aes(color = clusters), size = 0.7) +
  scale_color_manual(values = custom_colors) +
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
  ggtitle("UMAP by Cluster") +
  geom_text(data = centroids, aes(x = UMAP_1, y = UMAP_2, label = clusters), color = "black", size = 5, fontface = "bold")

# Legend ----
data <- data.frame(
  y = length(clusterNum):1,  # Posiciones y para los cuadros
  color = clusterCols,
  number = clusterNum
)

# Crear el plot de la leyenda de clusters
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
legendClusters

# Combining all plots----
comb.1 <- wrap_plots(umap_clusters, ab_sample_clusters, legendClusters)
comb.1

# Pathway ploting ####
gsea_res <- list.files("results/GSEA_scanRes/", pattern = "\\.csv$", full.names = TRUE) %>%
  setNames(gsub("2026-02-03_SSc_CD4_GSEA_|\\.csv$", "", basename(.))) %>%
  lapply(read_csv)
g_files <- list.files(file.path("results/GSEA_scanRes/"), pattern = "*.csv", full.names = TRUE)

clean_gsea_name <- function(filename, prefix) {
  basename <- tools::file_path_sans_ext(basename(filename))
  cleaned <- sub("^2026-02-03_SSc_CD4_GSEA_", "", basename)
}
names(gsea_res) <- vapply(g_files, clean_gsea_name, character(1))
gsea_res_df <- do.call(rbind, gsea_res)

cluster_names <- paste0("cluster_",c(0:11))
## Pathway selection ####

hallmark_pseub <- read.csv("ShinyCD4/data/GSEA_TFactivity/2026-01-31_SSc_CD4_GSEA_pseudobulk_Hallmark.csv")
reactome_pseub <- read.csv("ShinyCD4/data/GSEA_TFactivity/2026-01-31_SSc_CD4_GSEA_pseudobulk_C2.csv")

hallmark_pseub_pos <- hallmark_pseub %>%
  filter(padj < 0.01, NES > 0) %>%
  arrange(NES) %>%
  slice_head(n=6)
hallmark_pseub_neg <- hallmark_pseub %>%
  filter(padj < 0.01, NES < 0) %>%
  arrange(NES) %>%
  slice_head(n=6)
hallmark_pseub_fil <- rbind(hallmark_pseub_pos,hallmark_pseub_neg)


reactome_pseub_pos <- reactome_pseub %>%
  filter(padj < 0.01, NES > 0) %>%
  arrange(NES) %>%
  slice_head(n=6)
reactome_pseub_neg <- reactome_pseub %>%
  filter(padj < 0.01, NES < 0) %>%
  arrange(NES) %>%
  slice_head(n=6)
reactome_pseub_fil <- rbind(reactome_pseub_pos,reactome_pseub_neg)

pseudo_path <- rbind(hallmark_pseub_fil,reactome_pseub_fil)
pseudo_path_v <- as.vector(pseudo_path$pathway)

#
## Reactome ####
gsea_reactome <- gsea_res_df %>%
  filter(pathway %in% pseudo_path_v) %>%
  select(-...1) %>% 
  filter(grepl("^REACTOME", pathway)) %>%
  filter(padj < 0.1)

gsea_reactome_wide <- gsea_reactome %>%
  select(pathway,NES,cluster) %>%
  pivot_wider(
    names_from = cluster,
    values_from = NES
  )

gsea_reactome_wide <- gsea_reactome_wide[, c("pathway",0:11)]
colnames(gsea_reactome_wide) <- c("pathway",cluster_names)

top_pathways <- gsea_reactome_wide %>%
  ungroup() %>%
  mutate(n_na = rowSums(is.na(across(where(is.numeric))))) %>%
  arrange(n_na) #%>%
#slice_head(n = 20)



top_path_longer <- top_pathways %>%
  pivot_longer(
    cols = cluster_names,
    names_to = "cluster",
    values_to = "NES"
  ) %>%
  mutate(
    pathway = factor(
      pathway,
      levels = top_pathways %>%
        arrange(desc(n_na)) %>%
        pull(pathway)
    )
  )



ggplot(top_path_longer, aes(x = factor(cluster,levels = cluster_names), y = pathway)) +
  geom_point(aes(size = abs(NES), fill = ifelse(NES > 0, "NES > 0", "NES < 0")), 
             shape = 21, stroke = 0.5) +
  theme_bw() + 
  theme(
    axis.title.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right",
    plot.title = element_text(hjust = .65),
    axis.text.x = element_text(size = 9, face = "bold",angle = 45, hjust = 1),
    axis.text.y = element_text(size = 9)
  )+
  labs(x="", y = "",title = "Reactome Pathways across Clusters", size = "abs(NES)",fill="Direction")+
  scale_fill_manual(values = c("NES > 0" = "#F77663", "NES < 0" = "#3D799E")) +
  scale_size_continuous(name = "NES - FDR 10%")


## Hallmark ####
#group_by(pathway,NES,cluster)
#c("0","1","2","3","4","5","6","7","8","9","10","11")

gsea_hallmark <- gsea_res_df %>%
  filter(pathway %in% pseudo_path_v) %>%
  select(-...1) %>% 
  filter(grepl("^HALLMARK", pathway)) %>%
  filter(padj < 0.1)

gsea_hallmark_wide <- gsea_hallmark %>%
  select(pathway,NES,cluster) %>%
  pivot_wider(
    names_from = cluster,
    values_from = NES
  )

gsea_hallmark_wide <- gsea_hallmark_wide[, c("pathway",0:11)]
colnames(gsea_hallmark_wide) <- c("pathway",cluster_names)

top_pathways_h <- gsea_hallmark_wide %>%
  ungroup() %>%
  mutate(n_na = rowSums(is.na(across(where(is.numeric))))) %>%
  arrange(n_na) %>%
  slice_head(n = 20)

top_pathways_h$pathway <- sub("^HALLMARK_", "", top_pathways_h$pathway)


top_path_longer_h <- top_pathways_h %>%
  pivot_longer(
    cols = cluster_names,
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


ggplot(top_path_longer_h, aes(x = factor(cluster,levels = cluster_names), y = pathway)) +
  geom_point(aes(size = abs(NES), fill = ifelse(NES > 0, "NES > 0", "NES < 0")), 
             shape = 21, stroke = 0.5) +
  theme_bw() + 
  theme(
    axis.title.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "top",
    plot.title = element_text(hjust = .65),
    axis.text.x = element_text(size = 9, face = "bold",angle = 45, hjust = 1),
    axis.text.y = element_text(size = 9)
  )+
  labs(x="", 
       y = "",
       title = "Hallmark Pathways across Clusters", 
       size = "abs(NES)",
       fill="Direction")+
  scale_fill_manual(values = c("NES > 0" = "#F77663", 
                               "NES < 0" = "#3D799E")) +
  scale_size_continuous(name = "NES - FDR 10%")

# final plot
ggplot(
  top_path_longer_h %>% filter(!is.na(NES)),
  aes(x = factor(cluster, levels = cluster_names), y = pathway)
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
    plot.title = element_text(hjust = 1.6,size = 14),
    axis.text.x = element_text(size = 9, face = "bold", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 9),
    legend.position = "top",
    legend.justification = c(1.25, 0),
    legend.box = "vertical",
    legend.spacing.y = unit(0.1, "cm"),
    legend.box.spacing = unit(0.1, "cm"),
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0)
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

#
# version 2
ggplot(
  top_path_longer_h %>% filter(!is.na(NES)),
  aes(x = factor(cluster, levels = cluster_names), y = pathway)
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
    axis.text.x = element_text(size = 9, face = "bold", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 9)
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
#

# DEA quantification ####
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

# 1. Ruta principal
base_path <- "ShinyCD4/data/dea_sc/"

# 3. Leer subcarpeta de cluster
cluster_files <- list.files(file.path(base_path, "cluster"), pattern = "*.csv", full.names = TRUE)
cluster_list <- lapply(cluster_files, read.csv)
names(cluster_list) <- vapply(cluster_files, clean_dea_name, character(1))

cluster_id <- names(cluster_list)

dea_clusters <- lapply(names(cluster_list), function(cl) {
  scanpy_to_dea(cluster_list[[cl]],cl)
})

names(dea_clusters) <- names(cluster_list)

dea_cluster_df <- do.call(rbind,dea_clusters)
cluster_order <- paste0("Cluster",c(0:11))

dea_sig <- dea_cluster_df %>%
  filter(padj < 0.05) %>%        # o < 0.1 si quieres
  mutate(direction = ifelse(logFC > 0, "Up", "Down"))

dea_counts <- dea_sig %>%
  count(cluster, direction) %>%
  mutate(
    n_signed = ifelse(direction == "Up", n, -n)
  )

dea_counts$cluster <- factor(dea_counts$cluster, levels = rev(cluster_order))

ggplot(dea_counts, aes(x = cluster, y = n_signed, fill = direction)) +
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
    legend.position = "right"
  ) +
  labs(
    y = "Number of significant genes",
    fill = "Direction",
    title = "Significant DE genes per cluster in SSc"
  )


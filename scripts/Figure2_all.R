#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(pheatmap)
  library(readr)
  library(stringr)
  library(patchwork)
  library(fgsea)
  library(msigdbr)
  library(cowplot)
  library(ggh4x)
  library(Matrix)
})

# ------------------------------------------------
# Load metadata
# ------------------------------------------------

meta <- read.csv("data/20260124_SSc_CD4_meta_v3.csv", row.names = 1)
meta$clusters <- as.character(meta$R09)

cluster_order  <- 0:11
cluster_names  <- paste0("cluster_", cluster_order)

# ------------------------------------------------
# Color palettes
# ------------------------------------------------

cluster_colors <- c(
  "0"="#CDB4DB","1"="#F78181","2"="#82CFFD","3"="#8FCFD1",
  "4"="#AEB6BF","5"="#F7ACB4","6"="#FFB347","7"="#7FFFD4",
  "8"="#F4D03F","9"="#F8D5E1","10"="#FFA07A","11"="#2ECC71"
)

condition_colors <- c(
  "CTRL"="#D8BFD8",
  "SSc"="#C3E5C3"
)

direction_colors <- c(
  "Up"="#F77663",
  "Down"="#3D799E"
)

# ------------------------------------------------
# Abundance by sample
# ------------------------------------------------

abundance_sample <- meta %>%
  count(Sample, condition, R09) %>%
  group_by(Sample) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# ------------------------------------------------
# Abundance by condition
# ------------------------------------------------

abundance_condition <- meta %>%
  count(condition, R09) %>%
  group_by(condition) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()


ggplot(
  abundance_sample,
  aes(x = factor(condition, levels = c("CTRL","SSc")),
      y = prop,
      fill = condition)
) +
  geom_boxplot(width = 0.5, alpha = 0.75) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.6) +
  facet_wrap2(
    ~ factor(R09, levels = cluster_order),
    nrow = 4,
    strip = strip_themed(
      background_x = elem_list_rect(
        fill = cluster_colors,
        colour = "black"
      )
    )
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(size = 10),
    axis.title = element_blank()
  ) +
  scale_fill_manual(values = condition_colors)


ggplot(
  abundance_condition,
  aes(fill = factor(R09, levels = cluster_order),
      x = condition,
      y = prop)
) +
  geom_bar(stat = "identity", color = "black") +
  theme_classic() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_fill_manual(values = cluster_colors)


# ------------------------------------------------
# UMAP with cluster labels
# ------------------------------------------------

centroids <- meta %>%
  group_by(clusters) %>%
  summarize(
    UMAP_1 = mean(UMAP_1),
    UMAP_2 = mean(UMAP_2)
  )

ggplot(meta, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = clusters), size = 0.7) +
  scale_color_manual(values = cluster_colors) +
  theme_bw() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none",
    aspect.ratio = 1
  ) +
  geom_text(
    data = centroids,
    aes(label = clusters),
    size = 5,
    fontface = "bold"
  )

# ------------------------------------------------
# Differential expression analysis results
# ------------------------------------------------

scanpy_to_dea <- function(df, cluster_id) {
  out <- data.frame(
    gene    = df$SSc_names,
    logFC   = df$SSc_logfoldchanges,
    pval    = df$SSc_pvals,
    padj    = df$SSc_pvals_adj,
    zscore  = df$SSc_scores,
    cluster = cluster_id,
    stringsAsFactors = FALSE
  )
  out[!is.na(out$gene), ]
}

base_path <- "ShinyCD4/data/dea_sc/cluster"
cluster_files <- list.files(base_path, pattern="*.csv", full.names=TRUE)

cluster_list <- lapply(cluster_files, read.csv)
names(cluster_list) <- tools::file_path_sans_ext(basename(cluster_files))

dea_clusters <- lapply(names(cluster_list), function(cl) {
  scanpy_to_dea(cluster_list[[cl]], cl)
})

dea_cluster_df <- bind_rows(dea_clusters)

dea_counts <- dea_cluster_df %>%
  filter(padj < 0.05) %>%
  mutate(direction = ifelse(logFC > 0, "Up", "Down")) %>%
  count(cluster, direction) %>%
  mutate(n_signed = ifelse(direction == "Up", n, -n))

dea_counts$cluster <- factor(
  dea_counts$cluster,
  levels = rev(paste0("Cluster", cluster_order))
)

ggplot(dea_counts, aes(x = cluster, y = n_signed, fill = direction)) +
  geom_col(color = "black") +
  scale_fill_manual(values = direction_colors) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(
    y = "Number of significant genes",
    fill = "Direction",
    title = "Significant DE genes per cluster in SSc"
  )

# ------------------------------------------------
# GSEA results
# ------------------------------------------------

gsea_files <- list.files("results/GSEA_scanRes/",
                         pattern="\\.csv$",
                         full.names=TRUE)

gsea_list <- lapply(gsea_files, read_csv)
names(gsea_list) <- tools::file_path_sans_ext(basename(gsea_files))

gsea_res_df <- bind_rows(gsea_list)

# Select pathways based on pseudobulk filtering
hallmark_pseudobulk <- read.csv(
  "ShinyCD4/data/GSEA_TFactivity/2026-01-31_SSc_CD4_GSEA_pseudobulk_Hallmark.csv"
)

reactome_pseudobulk <- read.csv(
  "ShinyCD4/data/GSEA_TFactivity/2026-01-31_SSc_CD4_GSEA_pseudobulk_C2.csv"
)

pseudo_pathways <- bind_rows(
  hallmark_pseudobulk %>% filter(padj < 0.01) %>% slice_head(n=6),
  reactome_pseudobulk %>% filter(padj < 0.01) %>% slice_head(n=6)
) %>% pull(pathway)

# ------------------------------------------------
# IFN signature scores Heatmap
# ------------------------------------------------

IFNa_present <- IFNa_genes[IFNa_genes %in% rownames(exprs_norm)]
IFNg_present <- IFNg_genes[IFNg_genes %in% rownames(exprs_norm)]

meta$IFNa_score <- Matrix::colMeans(exprs_norm[IFNa_present, , drop=FALSE])
meta$IFNg_score <- Matrix::colMeans(exprs_norm[IFNg_present, , drop=FALSE])

score_df <- meta %>%
  group_by(clusters, condition) %>%
  summarise(
    IFNa = mean(IFNa_score),
    IFNg = mean(IFNg_score),
    .groups="drop"
  )

score_wide <- score_df %>%
  pivot_wider(
    names_from = condition,
    values_from = c(IFNa, IFNg),
    names_sep="_"
  )

score_wide$clusters <- factor(score_wide$clusters,
                              levels=as.character(cluster_order))

score_wide <- arrange(score_wide, clusters)

score_mat <- as.matrix(score_wide[, -1])
rownames(score_mat) <- paste0("cluster_", score_wide$clusters)

score_mat <- score_mat[, c(
  "IFNa_CTRL","IFNa_SSc",
  "IFNg_CTRL","IFNg_SSc"
)]

annotation_col <- data.frame(
  Signature=c("IFNa","IFNa","IFNg","IFNg"),
  Condition=c("CTRL","SSc","CTRL","SSc")
)

rownames(annotation_col) <- colnames(score_mat)

annotation_colors <- list(
  Signature=c(IFNa="#5E548E", IFNg="#1B998B"),
  Condition=condition_colors
)

gaps_row <- 2

ph <- pheatmap(
  t(score_mat),
  color=viridis(100),
  cluster_rows=FALSE,
  cluster_cols=FALSE,
  annotation_row=annotation_col,
  annotation_colors=annotation_colors,
  gaps_row=gaps_row,
  legend=FALSE,
  main="Cluster-level IFN-α and IFN-γ signature score"
)

min_val <- min(score_mat)
max_val <- max(score_mat)

legend_df <- data.frame(
  x=1,
  y=seq(min_val, max_val, length.out=100)
)

legend_plot <- ggplot(legend_df, aes(x=x, y=y, fill=y)) +
  geom_tile() +
  scale_fill_viridis_c(
    name="Mean of expression\nin group",
    limits=c(min_val, max_val),
    guide=guide_colorbar(
      direction="horizontal",
      title.position="top",
      barwidth=unit(5,"cm"),
      barheight=unit(0.4,"cm")
    )
  ) +
  theme_void() +
  theme(
    legend.position="bottom",
    legend.title=element_text(size=10),
    legend.text=element_text(size=9)
  )

legend_only <- cowplot::get_legend(legend_plot)

final_plot <- cowplot::plot_grid(
  ph$gtable,
  legend_only,
  ncol=1,
  rel_heights=c(4,1)
)

final_plot
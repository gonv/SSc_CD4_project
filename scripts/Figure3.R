#!/usr/bin/env Rscript

# ------------------------------------------------------------
# Figure 3 generation script
# CD4 T cell single-cell RNA-seq analysis in SSc
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(pheatmap)
  library(viridis)
  library(patchwork)
  library(cowplot)
  library(ggh4x)
  library(ggrepel)
  library(Matrix)
})

# ------------------------------------------------------------
# Create output directory
# ------------------------------------------------------------

dir.create("figures/plots4figures/Figure3", recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Load metadata
# ------------------------------------------------------------

meta <- read.csv("data/20260124_SSc_CD4_meta_v3.csv",
                 row.names = 1,
                 stringsAsFactors = FALSE)

metath <- read.csv("data/SSc_CD4_annotTh_metadata_with_umap.csv",
                   row.names = 1,
                   stringsAsFactors = FALSE)

metatreg <- read.csv("data/SSc_CD4_annotTreg_metadata_with_umap.csv",
                     row.names = 1,
                     stringsAsFactors = FALSE)

# Remove V gene and chain columns from Treg metadata
IR_V.cols <- grep("^IR_V", colnames(metatreg), value = TRUE)
chain.cols <- grep("chain", colnames(metatreg), value = TRUE)

metatreg <- metatreg %>%
  select(-all_of(IR_V.cols)) %>%
  select(-all_of(chain.cols))

# ------------------------------------------------------------
# Define annotations and color palette
# ------------------------------------------------------------

cellt <- c("Tn","Tcm","Tem","Treg","Temra")

annot_layer2 <- c("Th1","Th2","Th17","Th22","Thf")

annot_treg <- c("nTreg","iTreg","eTreg")

custom_colors <- c(
  "Tn"= "#a9dfbb",
  "Tcm"= "#6cd0c4",
  "Tem"= "#3faccd",
  "Treg"= "#fc6daf",
  "Temra"= "#ff0000",
  "Th17"= "#d7ff55",
  "Th2"= "#FF9F55",
  "Thf"= "#559FFF",
  "Th1"= "#FFC255",
  "Th22"= "#55FF9F",
  "nTreg" = "#6A4C93",
  "iTreg" = "#008B8B",
  "eTreg" = "#B22222",
  "CTRL" = "#D8BFD8",
  "SSc" = "#C3E5C3"
)

# ------------------------------------------------------------
# Layer 1 - Cell type abundance
# ------------------------------------------------------------

abundance_sample <- meta %>%
  mutate(
    cell_type = case_when(
      cell_type == "CTL" ~ "Temra",
      TRUE ~ cell_type
    )
  ) %>%
  count(Sample, condition, cell_type) %>%
  group_by(Sample) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

abundance_condition <- meta %>%
  mutate(
    cell_type = case_when(
      cell_type == "CTL" ~ "Temra",
      TRUE ~ cell_type
    )
  ) %>%
  count(condition, cell_type) %>%
  group_by(condition) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# ------------------------------------------------------------
# Boxplot (Figure 3C)
# ------------------------------------------------------------

fig3c <- ggplot(abundance_sample, aes(x = factor(condition, levels=c("CTRL","SSc")), y = prop, fill = condition)) +
  geom_boxplot(width = 0.5, alpha = 0.75) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.6) +
  facet_wrap2(
    ~ factor(cell_type, levels = cellt),
    scales = "free",
    nrow = 1,
    strip = strip_themed(
      background_x = elem_list_rect(
        fill = c( "#a9dfbb", "#6cd0c4", "#3faccd", "#fc6daf","#ff0000" ),
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
    strip.text = element_text( face = "bold", size = 12)
  ) +
  scale_fill_manual(
    values = c("#D8BFD8", "#C3E5C3")
  )

ggsave("figures/plots4figures/Figure3/figure3C.png",
       fig3c,
       width = 9.2,
       height = 1.5,
       dpi = 300)

# ------------------------------------------------------------
# UMAP annotation + Bar plot(Figure 3B)
# ------------------------------------------------------------

## UMAP
centroids <- meta %>%
  group_by(cell_type) %>%
  summarize(
    UMAP_1 = mean(UMAP_1),
    UMAP_2 = mean(UMAP_2),
    .groups = "drop"
  )

umap_annotation <- ggplot(meta, aes(x = UMAP_1, y = UMAP_2, color = cell_type)) +
  geom_point(color = "black", size = 0.9) +
  geom_point(aes(color = cell_type), size = 0.7) +
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
  ggtitle("T CD4+ cells annotation") +
  geom_text(data = centroids, aes(x = UMAP_1, y = UMAP_2, label = cell_type), color = "black", size = 4, fontface = "bold")

## Bar plots
ab_sample_annot <- ggplot(abundance_sample, aes(fill=factor(cell_type, levels =rev(cellt)), x=Sample, y=prop)) +
  geom_bar(stat = "identity", color='black')+
  labs(fill="",
       x=NULL,
       y="Fraction of cells")+
  theme_void()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #legend.key = element_rect(fill = "lightblue", color = NA),
        legend.text = element_text(size = 11, face = "bold"),
        legend.key.size = unit(0.5, "cm"),
        legend.key.width = unit(0.5,"cm"), 
        axis.ticks = element_blank(),
        axis.text=element_text(size=11),
        axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1),
        axis.title=element_blank(), 
        plot.margin = margin(1,0.1,0.1,0.1, "cm")
  )+
  scale_fill_manual(guide = guide_legend(reverse=TRUE), 
                    values=custom_colors)

## Combine plots (Figure 3B)
comb.1 <- wrap_plots(umap_annotation, ab_sample_annot)

ggsave("figures/plots4figures/Figure3/figure3B.png",
       comb.1,
       width = 6.7,
       height = 3.5,
       dpi = 300)

# ------------------------------------------------------------
# Layer 2 - T helper annotation
# ------------------------------------------------------------
rm(abundance_sample,abundance_condition)

abundance_sample <- metath %>%
  count(Sample, condition, Th) %>%
  group_by(Sample) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

abundance_condition <- metath %>%
  dplyr::count(condition, Th) %>%
  group_by(condition) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()


# ------------------------------------------------------------
# Boxplot (Figure 3F)
# ------------------------------------------------------------

fig3F <- ggplot(abundance_sample, aes(x = factor(condition, levels=c("CTRL","SSc")), y = prop, fill = condition)) +
  geom_boxplot(width = 0.5, alpha = 0.75) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.6) +
  facet_wrap2(
    ~ factor(Th, levels = annot_layer2),
    scales = "free",
    nrow = 5,
    strip = strip_themed(
      background_x = elem_list_rect(
        fill = c( "#FFC255", "#FF9F55", "#d7ff55", "#55FF9F","#559FFF" ),
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
    strip.text = element_text( face = "bold", size = 12)
  ) +
  scale_fill_manual(
    values = c("#D8BFD8", "#C3E5C3")
  )

ggsave(plot = fig3F, 
       filename = "figures/plots4figures/Figure3/figure3F.png",
       width = 1.5,
       height = 6,
       dpi = 300)

# ------------------------------------------------------------
# UMAP T helper annotation (Figure 3E)
# ------------------------------------------------------------

centroids.1 <- metath %>%
  group_by(Th) %>%
  summarize(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2))

umap_sub_th <- ggplot(metath, aes(x = UMAP_1, y = UMAP_2, color = Th)) +
  geom_point(color = "black", size = 0.9) +
  geom_point(aes(color = Th), size = 0.7) +
  scale_color_manual(values = custom_colors) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, 
                              face = "bold", 
                              size = 14),
    axis.title = element_text(face = "bold", 
                              size = 12),
    axis.text = element_blank(),
    legend.text = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = "black", size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 14)
  ) +
  ggtitle("T helper annotation") +
  geom_text(data = centroids.1, aes(x = UMAP_1, y = UMAP_2, label = Th), color = "black", size = 5, fontface = "bold")

centroids.2 <- meta %>%
  group_by(Th) %>%
  summarize(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2))

umap_th <- ggplot(meta, aes(x = UMAP_1, y = UMAP_2, color = Th)) +
  geom_point(size = 1) +
  scale_color_manual(values = custom_colors) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, 
                              face = "bold", 
                              size = 12),
    axis.title = element_text(face = "bold", 
                              size = 12),
    axis.text = element_blank(),
    legend.text = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = "black", size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 14)
  ) +
  ggtitle("") +
  geom_text_repel(
    data = centroids.2,
    aes(x = UMAP_1, y = UMAP_2, label = Th),
    size = 6,
    fontface = "bold",
    color = "black",
    box.padding = 0.5,
    point.padding = 0.5,
    segment.color = "black",
    max.overlaps = Inf
  )

fig3E <- wrap_plots(umap_sub_th,umap_th, ncol = 1)

ggsave(plot = fig3E, 
       filename = "figures/plots4figures/Figure3/figure3E.png",
       width = 2.9,
       height = 5.5,
       dpi = 300)

# ------------------------------------------------------------
# Heatmap (Figure 3G)
# ------------------------------------------------------------

exprs_norm <- readRDS("ShinyCD4/data/20260124_SSc_CD4_counts_norm.rds")

genes_th <- c(
  "VIM","ANXA1","SELL","CCR4","CD48","KLRB1",
  "CCL20","LTK","ABCB1","CTSH","PTPN13",
  "CCL5","HLA-DRB1","CD74","HLA-DPB1"
)

cells_tcm <- meta$cell_type == "Tcm"

meta_th <- meta[cells_tcm,]
exprs_tcm <- exprs_norm[,cells_tcm]

genes_present <- genes_th[genes_th %in% rownames(exprs_tcm)]

exprs_sub <- exprs_tcm[genes_present,,drop=FALSE]

meta_th$Th_condition <- paste(meta_th$Th,
                              meta_th$condition,
                              sep="_")

groups <- factor(meta_th$Th_condition)

design <- sparse.model.matrix(~0+groups)

colnames(design) <- levels(groups)

exprs_sum <- exprs_sub %*% design

group_sizes <- colSums(design)

exprs_mean <- t(t(exprs_sum) / group_sizes)

exprs_mean <- as.matrix(exprs_mean)

# Row scaling
scale_rows <- function(mat){
  m <- rowMeans(mat)
  s <- apply(mat,1,sd)
  s[s==0] <- 1
  sweep(sweep(mat,1,m,"-"),1,s,"/")
}

exprs_scaled <- scale_rows(exprs_mean)

exprs_scaled[exprs_scaled > 2] <- 2
exprs_scaled[exprs_scaled < -2] <- -2

heat_th <- pheatmap(
  exprs_scaled,
  color = viridis(100),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  border_color = "black",
  angle_col = 90,
  fontsize_col = 10, 
  filename = "figures/plots4figures/Figure3/figure3G.png", 
  width = 3, 
  height = 5.2
)

# ------------------------------------------------------------
# Layer 3 - T regulatory cells annotation
# ------------------------------------------------------------
rm(abundance_sample,abundance_condition)

abundance_sample <- metatreg %>%
  count(Sample, condition, ct_Treg) %>%
  group_by(Sample) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

abundance_condition <- metatreg %>%
  dplyr::count(condition, ct_Treg) %>%
  group_by(condition) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# ------------------------------------------------------------
# UMAP T regulatory cells annotation (Figure 3J)
# ------------------------------------------------------------

centroids.3 <- metatreg %>%
  group_by(ct_Treg) %>%
  summarize(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2))

umap_treg_sub <- ggplot(metatreg, 
                        aes(x = UMAP_1, y = UMAP_2, color = ct_Treg)
                        ) +
  geom_point(color = "black", size = 0.9) +
  geom_point(aes(color = ct_Treg), size = 0.7) +
  scale_color_manual(values = custom_colors) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_blank(),
    legend.text = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = "black", size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 14)
  ) +
  ggtitle("T regulatory cells\nannotation") +
  geom_text(data = centroids.3, 
            aes(x = UMAP_1, y = UMAP_2, label = ct_Treg), 
            color = "white", 
            size = 6, 
            fontface = "bold")

## Merging metas
metatreg$cell <- rownames(metatreg)
meta_temp <- metatreg %>% select(cell, ct_Treg)

meta_join <- full_join(meta, meta_temp, by = "cell")
meta_join <- meta_join %>% 
  mutate(ct_Treg = coalesce(ct_Treg, cell_type))

centroids.4 <- meta_join %>%
  group_by(ct_Treg) %>%
  summarize(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2))

umap_treg <- ggplot(meta_join, aes(x = UMAP_1, y = UMAP_2, color = ct_Treg)) +
  geom_point(size = 1) +
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
    strip.text = element_text(face = "bold", size = 14)
  ) +
  ggtitle("") +
  geom_text_repel(
    data = centroids.4,
    aes(x = UMAP_1, y = UMAP_2, label = ct_Treg),
    size = 6,
    fontface = "bold",
    color = "black",
    box.padding = 0.5,
    point.padding = 0.5,
    segment.color = "black",
    max.overlaps = Inf
  )

fig3I <- wrap_plots(umap_treg_sub,umap_treg, ncol = 1)

ggsave(plot = fig3I, 
       filename = "figures/plots4figures/Figure3/figure3I.png",
       width = 2.2,
       height = 5,
       dpi = 300)

# ------------------------------------------------------------
# Boxplot (Figure 3J)
# ------------------------------------------------------------

fig3J <- ggplot(abundance_sample, 
                aes(x = factor(condition, levels=c("CTRL","SSc")), y = prop, fill = condition)
                ) +
  geom_boxplot(width = 0.5, alpha = 0.75) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.6) +
  facet_wrap2(
    ~ factor(ct_Treg, levels = annot_treg),
    scales = "free",
    nrow = 1,
    strip = strip_themed(
      background_x = elem_list_rect(
        fill = c( "#6A4C93", "#008B8B", "#B22222" ),
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
    strip.text = element_text(color = "white", face = "bold", size = 12)
  ) +
  scale_fill_manual(
    values = c("#D8BFD8", "#C3E5C3")
  )

ggsave(plot = fig3J, 
       filename = "figures/plots4figures/Figure3/figure3J.png",
       width = 4.4,
       height = 1,
       dpi = 300)

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

# Colors and vectors----
cellt <-c("Tn","Tcm","Tem","Treg","Temra")
annot_layer2 <- c("Th1","Th2","Th17","Th22","Thf")
annot_treg <- c("nTreg","iTreg","eTreg")

custom_colors <- c(
  "0" = "#CDB4DB", "1" = "#F78181", "2" = "#82CFFD", "6" = "#FFB347", "4" = "#AEB6BF",
  "5" = "#F7ACB4", "7" = "#7FFFD4", "10" = "#FFA07A", "8" = "#F4D03F", "9" = "#F8D5E1",
  "3" = "#8FCFD1", "11" = "#2ECC71", "CTRL" = "#D8BFD8", "SSc" = "#C3E5C3",
  "Tn"= "#a9dfbb", "Tcm"= "#6cd0c4", "Treg"= "#fc6daf", "Tem"= "#3faccd",
  "Temra"= "#ff0000", "Th17"= '#d7ff55', "Th2"= '#FF9F55', "Thf"= '#559FFF',
  "Th1_GZMK"= '#55FF9F', "Th1_GZMB"= '#FFC255',"Th22"='#55FF9F',"Th1"='#FFC255',
  "nTreg" = "#6A4C93","iTreg" = "#008B8B","eTreg" ="#B22222"
)

# Layer 1 ----
## Abundance ----
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
  dplyr::count(condition, cell_type) %>%
  group_by(condition) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

## Plots ----
### Box plots ----
ggplot(abundance_sample, aes(x=factor(condition, levels=c("CTRL","SSc")), y=prop, fill=condition)) + 
  facet_wrap(~cell_type, ncol = 5)+
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
                    values=colorCT)

ggplot(abundance_sample, aes(x = factor(condition, levels=c("CTRL","SSc")), y = prop, fill = condition)) +
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
    plot.margin = unit(c(1, 1, 1, 1), "mm")
  ) +
  scale_fill_manual(
    values = c("#D8BFD8", "#C3E5C3")
  )
#
### Bar plots ----
# By condition
ggplot(abundance_condition, aes(fill=factor(cell_type, levels =cellt), x=condition, y=prop)) +
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
ggplot(abundance_sample, aes(fill=factor(cell_type, levels =rev(cellt)), x=Sample, y=prop)) +
  geom_bar(stat = "identity", color='black')+
  labs(fill="Cell_Type",x=NULL,y="Fraction of cells")+
  theme_void()+
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
#
# Layer 2  ----
## Abundance ----
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

## Plots ----
### Box plots ----
ggplot(abundance_sample, aes(x=factor(condition, levels=c("CTRL","SSc")), y=prop, fill=condition)) + 
  facet_wrap(~Th, ncol = 5)+
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
    plot.margin = unit(c(1, 1, 1, 1), "mm")
  ) +
  scale_fill_manual(
    values = c("#D8BFD8", "#C3E5C3")
  )
#
### Bar plots ----
# By condition
ggplot(abundance_condition, aes(fill=factor(Th, levels =annot_layer2), x=condition, y=prop)) +
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
ggplot(abundance_sample, aes(fill=factor(Th, levels =rev(annot_layer2)), x=Sample, y=prop)) +
  geom_bar(stat = "identity", color='black')+
  labs(fill="Cell_Type",x=NULL,y="Fraction of cells")+
  theme_void()+
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
### UMAPs ----

centroids <- metath %>%
  group_by(Th) %>%
  summarize(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2))

umap_sub_th <- ggplot(metath, aes(x = UMAP_1, y = UMAP_2, color = Th)) +
  geom_point(color = "black", size = 0.9) +
  geom_point(aes(color = Th), size = 0.7) +
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
  ggtitle("UMAP by T helper annotation \n subclustering") +
  geom_text(data = centroids, aes(x = UMAP_1, y = UMAP_2, label = Th), color = "black", size = 5, fontface = "bold")

centroids <- meta %>%
  group_by(Th) %>%
  summarize(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2))

umap_th <- ggplot(meta, aes(x = UMAP_1, y = UMAP_2, color = Th)) +
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
  geom_text(data = centroids, aes(x = UMAP_1, y = UMAP_2, label = Th), color = "black", size = 5, fontface = "bold")

wrap_plots(umap_sub_th,umap_th, ncol = 1)
#
# Layer 3----
## Abundance ----
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

## Plots ----
### Box plots ----
ggplot(abundance_sample, aes(x=factor(condition, levels=c("CTRL","SSc")), y=prop, fill=condition)) + 
  facet_wrap(~ct_Treg, ncol = 5)+
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
    plot.margin = unit(c(1, 1, 1, 1), "mm")
  ) +
  scale_fill_manual(
    values = c("#D8BFD8", "#C3E5C3")
  )
#
### Bar plots ----
# By condition
ggplot(abundance_condition, aes(fill=factor(ct_Treg, levels =annot_treg), x=condition, y=prop)) +
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
ggplot(abundance_sample, aes(fill=factor(ct_Treg, levels =rev(annot_treg)), x=Sample, y=prop)) +
  geom_bar(stat = "identity", color='black')+
  labs(fill="Cell_Type",x=NULL,y="Fraction of cells")+
  theme_void()+
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
### UMAPs ----
centroids <- metatreg %>%
  group_by(ct_Treg) %>%
  summarize(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2))

ggplot(metatreg, aes(x = UMAP_1, y = UMAP_2, color = Th)) +
  geom_point(color = "black", size = 0.9) +
  geom_point(aes(color = ct_Treg), size = 0.7) +
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
  geom_text(data = centroids, aes(x = UMAP_1, y = UMAP_2, label = ct_Treg), color = "black", size = 5, fontface = "bold")



centroids <- metatreg %>%
  group_by(ct_Treg) %>%
  summarize(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2))

umap_treg_sub <- ggplot(metatreg, aes(x = UMAP_1, y = UMAP_2, color = ct_Treg)) +
  geom_point(color = "black", size = 0.9) +
  geom_point(aes(color = ct_Treg), size = 0.7) +
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
  ggtitle("UMAP by Treg annotation \n subclustering") +
  geom_text(data = centroids, aes(x = UMAP_1, y = UMAP_2, label = ct_Treg), color = "black", size = 5, fontface = "bold")

### Merging metas ----
metatreg$cell <- rownames(metatreg)
meta_temp <- metatreg %>% select(cell, ct_Treg)

meta_join <- full_join(meta, meta_temp, by = "cell")
meta_join <- meta_join %>% 
  mutate(ct_Treg = coalesce(ct_Treg, cell_type))

centroids <- meta_join %>%
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
  geom_text(data = centroids, aes(x = UMAP_1, y = UMAP_2, label = ct_Treg), color = "black", size = 5, fontface = "bold")
umap_treg

wrap_plots(umap_treg_sub,umap_treg, ncol = 1)
#
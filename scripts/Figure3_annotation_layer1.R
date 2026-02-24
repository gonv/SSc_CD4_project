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

custom_colors <- c(
  "0" = "#CDB4DB", "1" = "#F78181", "2" = "#82CFFD", "6" = "#FFB347", "4" = "#AEB6BF",
  "5" = "#F7ACB4", "7" = "#7FFFD4", "10" = "#FFA07A", "8" = "#F4D03F", "9" = "#F8D5E1",
  "3" = "#8FCFD1", "11" = "#2ECC71", "CTRL" = "#D8BFD8", "SSc" = "#C3E5C3",
  "Tn"= "#a9dfbb", "Tcm"= "#6cd0c4", "Treg"= "#fc6daf", "Tem"= "#3faccd",
  "Temra"= "#ff0000", "Th17"= '#d7ff55', "Th2"= '#FF9F55', "Thf"= '#559FFF',
  "Th1_GZMK"= '#55FF9F', "Th1_GZMB"= '#FFC255',"Th22"='#55FF9F',"Th1"='#FFC255'
)

# Abundance ----
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

# Plots ----
## Box plots ----
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
## Bar plots ----
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
##  ----
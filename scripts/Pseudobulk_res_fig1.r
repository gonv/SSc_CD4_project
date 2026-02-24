setwd("/lab-share/IM-Gutierrez-e2/Public/gvillanueva/00_CD4/analysis/")

.libPaths("/home/ch258010/miniconda3/envs/condaR441/lib/R/library")

data_lib <- c("devtools","tidyverse","forcats","purrr","readr","tibble","stringr","data.table","DT")
pathway_lib <- c("msigdbr","fgsea")
TF_lib <- c("dorothea","decoupleR","viper","OmnipathR")
RNAseq_lib <- c("edgeR","lme4")
Seurat_lib <- c("Seurat","harmony","sctransform","MAST","scater")
Single_cell <- c("Matrix","parallel","magrittr","umap","presto","ggtree","aricode","genesorteR")
ploting_lib <- c("pheatmap","RColorBrewer","ggridges","patchwork","cowplot","viridis","gridExtra","ggpubr","NMF","purrr","ggrepel")
#ATAC_lib<-c("Signac","BSgenome.Hsapiens.UCSC.hg38","GenomicRanges","GenomeInfoDb","AnnotationHub")

print("Data wrangling packages")
lapply(data_lib, function(x) suppressMessages(suppressWarnings(require(x, character.only = TRUE))))
print("Pathway analysis packages")
lapply(pathway_lib, function(x) suppressMessages(suppressWarnings(require(x, character.only = TRUE))))
print("Transcription factors activity analysis packages")
lapply(TF_lib, function(x) suppressMessages(suppressWarnings(require(x, character.only = TRUE))))  
print("RNA-seq data management")
lapply(RNAseq_lib, function(x) suppressMessages(suppressWarnings(require(x, character.only = TRUE))))      
print("Single cell data management")
lapply(Seurat_lib, function(x) suppressMessages(suppressWarnings(require(x, character.only = TRUE))))
lapply(Single_cell, function(x) suppressMessages(suppressWarnings(require(x, character.only = TRUE))))     
print("Plotting data packages")
lapply(ploting_lib, function(x) suppressMessages(suppressWarnings(require(x, character.only = TRUE))))
#print("ATAC-seq analysis packages")
#lapply(ATAC_lib, function(x) suppressMessages(suppressWarnings(require(x, character.only = TRUE))))
       
set.seed("42")

print("Qian functions")
## Functions
source("/lab-share/IM-Gutierrez-e2/Public/gvillanueva/00_CD4/scripts/someFunctions_SR.R")
source("/lab-share/IM-Gutierrez-e2/Public/Qian/B_Atlas/analysis_scripts/Some_functions.r")

fig.size <- function (height, width) {
    options(repr.plot.height = height, repr.plot.width = width)
}

res2_ann <- read.csv("../results/dea_pseudobulk_condition_v6.csv")

res2_ann %>% dim
res2_ann %>% str
res2_ann %>% DT::datatable()

res2_ann$adj.P.Val <- p.adjust(res2_ann$p_value, method = "fdr")
res2_ann$significance <- case_when(
  res2_ann$adj.P.Val < 0.05 & res2_ann$Estimate > 0  ~ "Upregulated",
  res2_ann$adj.P.Val < 0.05 & res2_ann$Estimate < 0 ~ "Downregulated",
  TRUE ~ "Not_significant"
)
res2_ann$delabel <- NA
res2_ann$delabel[res2_ann$significance != "Not_significant"] <- res2_ann$gene[res2_ann$significance != "Not_significant"]

fig.size(8,10)
fdr_threshold <- -log10(max(res2_ann$p_value[res2_ann$adj.P.Val <= 0.05], na.rm = TRUE))
volplot <- ggplot(res2_ann, aes(x = Estimate, y = -log10(p_value),label = delabel)) +
  geom_point(alpha = 0.9, 
             shape=21, 
             colour="black", 
             stroke=0.3, 
             size=2.5, 
             aes(fill=significance)) + 
  theme_minimal() +  
  labs(title = "", 
       x = "Estimate", 
       y = "-log10(p-value)",
       fill="",
       color="") +
  geom_label_repel(data = filter(res2_ann, !is.na(delabel)),, color = "black",
                   max.overlaps = 10, force = 5, size = 3, show.legend = FALSE) +
  scale_fill_manual(values = c("Upregulated" = "#F77663", 
                               "Downregulated" = "#3D799E", 
                               "Not_significant" = "#f0f0f0")) +
  geom_hline(yintercept = fdr_threshold, linetype = "dashed", color = "blue")+
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")
volplot

res2_ann %>% DT::datatable()

res2_ann_summ <- res2_ann %>%
  select(significance) %>%
  group_by(significance) %>%
  summarise(n = n(), .groups = "drop")

res2_ann_summ

ls()

DE_pos <- res2_ann %>%
  filter(adj.P.Val < 0.05) %>%
  filter(Estimate > 0) %>%
  arrange(desc(Estimate)) %>%
  slice_head(n = 20)

DE_neg <- res2_ann %>%
  filter(adj.P.Val < 0.05) %>%
  filter(Estimate < 0) %>%
  arrange(Estimate) %>%
  slice_head(n = 20)

DE_combined <- bind_rows(DE_pos, DE_neg)

genes <- DE_combined$gene

res2_ann$fig <- ifelse(res2_ann$gene %in% genes, "TopGenes", "Other")
res2_ann$delabel <- NA
res2_ann$delabel[res2_ann$fig != "Other"] <- res2_ann$gene[res2_ann$fig != "Other"]

fig.size(8,10)
fdr_threshold <- -log10(max(res2_ann$p_value[res2_ann$adj.P.Val <= 0.05], na.rm = TRUE))

volplot <- ggplot(res2_ann, aes(x = Estimate, y = -log10(p_value),label = delabel)) +

  geom_point(alpha = 0.9, 
             shape=21, 
             colour="black", 
             stroke=0.3, 
             size=3, 
             aes(fill=significance)) + 

  theme_minimal() +
  theme(axis.title = element_text(size=16),
        legend.position = "bottom",
        legend.text = element_text(size=16),
        legend.title = element_blank())+

  labs(title = "", 
       x = "Estimate", 
       y = "-log10(p-value)",
       fill="",
       color="") +

  geom_label_repel(data = filter(res2_ann, !is.na(delabel)),, color = "black",
                   max.overlaps = 50, force = 5, size = 4, show.legend = FALSE) +

  scale_fill_manual(values = c("Upregulated" = "#F77663", 
                               "Downregulated" = "#3D799E", 
                               "Not_significant" = "#f0f0f0")) +

  geom_hline(yintercept = fdr_threshold, linetype = "dashed", color = "blue")+
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")

volplot

ggsave(volplot,filename = "../results/plots/20260203_SSc_CD4_PseudoBulk_volcano.png",width = 8,height = 8, dpi = 300)

temp_rec <- readRDS("2025-10-24_SSc_CD4_so_process_fnl_v3.rds")

temp_rec

counts_all <- temp_rec@assays$RNA@layers$counts

meta_ssc <- temp_rec@meta.data
umap_coor <- temp_rec@reductions$umap@cell.embeddings
pcs_coord <- temp_rec@reductions$pca@cell.embeddings
meta_ssc <- cbind(meta_ssc,pcs_coord[,1:2],umap_coor)
meta_ssc$cell <- rownames(meta_ssc)

new_names <- c("SSc_1","SSc_2","SSc_3","SSc_4","SSc_5","SSc_6","SSc_7","SSc_8","CTRL_1","CTRL_2","CTRL_3","CTRL_4","CTRL_5","CTRL_6","CTRL_7","CTRL_8")

sample_names <- data.frame(
    og = meta_ssc$Sample %>% unique,
    new = new_names
)

sample_names

counts_all %>% str

temp_rec@meta.data %>% head

meta_ssc$Sample %>% unique

# Generate pseudobulk matrix
#exprs_raw <- counts_cluster[, meta_ssc_cluster$cell]
exprs_norm <- counts_all %>% NormalizeDataSeurat() # for the heatmap

DE_pos <- res2_ann %>%
  filter(adj.P.Val < 0.05) %>%
  filter(Estimate > 0) %>%
  arrange(desc(Estimate)) %>%
  slice_head(n = 20)

DE_neg <- res2_ann %>%
  filter(adj.P.Val < 0.05) %>%
  filter(Estimate < 0) %>%
  arrange(Estimate) %>%
  slice_head(n = 20)

DE_combined <- bind_rows(DE_pos, DE_neg)

gene_plot <- unique(DE_combined$gene)
gene_plot %>% length

exp_heat <- exprs_norm[gene_plot,]
exp_heat <- as.data.frame(t(exp_heat))

exp_heat$sample <- as.character(meta_ssc$Sample)
dim(exp_heat)

exp_heat %>% head

exp_ave <- aggregate(exp_heat[, 1:(ncol(exp_heat)-1)], list(exp_heat$sample), mean)

colnames(exp_ave)[1] <- "sample"
exp_ave <- as.data.frame(t(exp_ave))
colnames(exp_ave) <- as.character(t(exp_ave[1,]))
exp_ave <- exp_ave[-1,]
row_names <- rownames(exp_ave)

exp_ave <- mutate_all(exp_ave, function(x) as.numeric(as.character(x)))
rownames(exp_ave) <- row_names
exp_ave <- as.matrix(exp_ave)
scale_rows <- function(x) t(scale(t(x)))
exp_ave_scale <- scale_rows(exp_ave) # Z-score
exp_ave_scale[exp_ave_scale > 2] <- 2
exp_ave_scale[exp_ave_scale < -2] <- -2
dim(exp_ave_scale)

exp_ave_scale %>% str
colnames(exp_ave_scale)
exp_ave_scale %>% head

col_info <- data.frame(
  sample = colnames(exp_ave_scale),
  condition = ifelse(grepl("^SSc", colnames(exp_ave_scale)), "SSc", "CTRL"),
  id = as.numeric(gsub(".*_", "", colnames(exp_ave_scale)))
)

col_info_ord <- col_info %>%
  dplyr::arrange(condition == "CTRL", id)

exp_ave_scale <- exp_ave_scale[, 
  colnames(exp_ave_scale)[
    order(
      grepl("^SSc", colnames(exp_ave_scale)),
      as.numeric(gsub(".*_", "", colnames(exp_ave_scale)))
    )
  ]
]

subset <- which(rownames(exp_ave_scale) %in% rownames(exp_ave_scale)) 
labels = rownames(exp_ave_scale)[subset]

row_annot <- data.frame(Link = labels)  # labels debe tener length = nrow(exp_ave_scale)
rownames(row_annot) <- rownames(exp_ave_scale)

exp_ave_scale %>% str

# Asignining new names to the heatmap
name_map <- setNames(sample_names$new, sample_names$og)

colnames(exp_ave_scale) <- ifelse(
  colnames(exp_ave_scale) %in% names(name_map),
  name_map[colnames(exp_ave_scale)],
  colnames(exp_ave_scale)
)

fig.size(10,9)
heatdea <- pheatmap(
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
  border_color = NA
)
heatdea

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
    filename = "../results/plots/20260203_SSc_CD4_PseudoBulk_heatmap.png",
    width = 9,
    height = 10, 
    dpi = 300
)

ggsave(volplot,filename = "../results/plots/20260203_SSc_CD4_PseudoBulk_volcano.png",width = 8,height = 8, dpi = 300)

c2_all <- msigdbr(species = "Homo sapiens", category = "C2")

# KEGG
kegg <- c2_all %>%
  filter(gs_subcat == "CP:KEGG") %>%
  split(.$gene_symbol, .$gs_name)

# Reactome
react <- c2_all %>%
  filter(gs_subcat == "CP:REACTOME") %>%
  split(.$gene_symbol, .$gs_name)

# WikiPathways
wp <- c2_all %>%
  filter(gs_subcat == "CP:WP") %>%
  split(.$gene_symbol, .$gs_name)

# Hallmark (similar a "canon" del script original)
hallmark <- msigdbr(species="human", category="H")
H_gene_sets <- msigdbr(species = "human", category = 'H')

H.ensembl.ls <- H_gene_sets %>% 
  select(gs_name, gene_symbol) %>% 
  group_by(gs_name) %>% 
  summarise(all.genes = list(unique(gene_symbol))) %>% 
  deframe()

# Ranking the genes for GSEA
ranking <- setNames(res2_ann$t_value, res2_ann$gene) # asign to each gene its OR

ranking %>% str

c2_gene_sets = msigdbr(species = "human", category = 'C2')

react <- c2_gene_sets %>%
  filter(gs_subcollection == "CP:REACTOME") %>%
  select(gs_name, gene_symbol) %>% 
  group_by(gs_name) %>% 
  summarise(all.genes = list(unique(gene_symbol))) %>% 
  deframe()

C2.ensembl.ls <- c2_gene_sets %>% 
  select(gs_name, gene_symbol) %>% 
  group_by(gs_name) %>% 
  summarise(all.genes = list(unique(gene_symbol))) %>% 
  deframe()

gsea.c2 <- fgseaMultilevel(pathways= react,
                           stats=ranking,  
                           minSize = 5, 
                           maxSize = 500)
head(gsea.c2[order(pval), ])

fgseaResTidy <- gsea.c2 %>%
  as_tibble() %>%
  arrange(desc(NES))

# Show in a nice table:
fgseaResTidy %>% 
  dplyr::select( -ES) %>% 
  arrange(padj) %>% 
  DT::datatable()

top_pos <- gsea.c2 %>%
  filter(padj < 0.1) %>%
  filter(NES > 0) %>%
  arrange(desc(NES)) %>%
  slice_head(n = 10)

top_neg <- gsea.c2 %>%
  filter(padj < 0.1) %>%
  filter(NES < 0) %>%
  arrange(NES) %>%
  slice_head(n = 10)

top_combined <- bind_rows(top_pos, top_neg)

fig.size(8,13)
gsea.c2.plot <- ggplot(top_combined, aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill = ifelse(NES > 0, "NES > 0", "NES < 0")), alpha = 0.8,width = 0.75,colour="black") +
    scale_fill_manual(values = c("NES > 0" = "#FF7F50", "NES < 0" = "#6495ED"), name = "FDR 10%") +
    coord_flip() +
    labs(x = "Pathway", y = "Normalized Enrichment Score",
         title = paste0("REACTOME top pathways in SSc vs CTRL")) + 
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12,face = "bold"),
      axis.title = element_text(size = 14),
      legend.position = "top",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )

gsea.c2.plot

fig.size(4,12)
plotEnrichment(react[["REACTOME_INTERFERON_ALPHA_BETA_SIGNALING"]], ranking) +
  labs(
    title = "REACTOME_INTERFERON_ALPHA_BETA_SIGNALING"
  )

fig.size(6,10)
plotEnrichment(react[["REACTOME_PYROPTOSIS"]], ranking) +
  labs(
    title = "Pyroptosis"
  )

fig.size(6,10)
plotEnrichment(react[["REACTOME_NGF_STIMULATED_TRANSCRIPTION"]], ranking) +
  labs(
    title = ""
  )

gsea.C2_out <- gsea.c2 %>%
    mutate(
        leadingEdge = purrr::map_chr(leadingEdge, ~ paste(.x, collapse = ";"))
          )

#write.csv(gsea.C2_out,"2026-01-31_SSc_CD4_GSEA_pseudobulk_C2.csv")

gsea.H <- fgseaMultilevel(pathways = H.ensembl.ls, 
                           stats=ranking,  
                           minSize = 10, 
                           maxSize = 500)

gsea.H_out <- gsea.H %>%
    mutate(
        leadingEdge = purrr::map_chr(leadingEdge, ~ paste(.x, collapse = ";"))
        )

#write.csv(gsea.H_out,"2026-01-31_SSc_CD4_GSEA_pseudobulk_Hallmark.csv")


fgseaResTidy <- gsea.H %>%
  as_tibble() %>%
  arrange(desc(NES))

# Show in a nice table:
fgseaResTidy %>% 
  dplyr::select(-ES) %>% 
  arrange(padj) %>% 
  DT::datatable()

fig.size(6,10)
plotEnrichment(H.ensembl.ls[["HALLMARK_INFLAMMATORY_RESPONSE"]], ranking) +
  labs(
    title = "Inflammatory response"
  )

plotEnrichment(H.ensembl.ls[["HALLMARK_TNFA_SIGNALING_VIA_NFKB"]], ranking) +
  labs(
    title = "TNFa via NKFB"
  )

fig.size(4,12)
plotEnrichment(H.ensembl.ls[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]], ranking) +
  labs(
    title = "HALLMARK_INTERFERON_GAMMA_RESPONSE"
  )

top_pos.H <- gsea.H %>%
  filter(padj < 0.1) %>%
  filter(NES > 0) %>%
  arrange(desc(NES)) %>%
  slice_head(n = 10)

top_neg.H <- gsea.H %>%
  filter(padj < 0.1) %>%
  filter(NES < 0) %>%
  arrange(NES) %>%
  slice_head(n = 10)

top_combined.H <- bind_rows(top_pos.H, top_neg.H)

fig.size(8,10)
gsea.H.plot <- ggplot(top_combined.H, aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill = ifelse(NES > 0, "NES > 0", "NES < 0")), alpha = 0.8,width = 0.75,colour="black") +
    scale_fill_manual(values = c("NES > 0" = "#FF7F50", "NES < 0" = "#6495ED"), name = "FDR 10%") +
    coord_flip() +
    labs(x = "Pathway", y = "Normalized Enrichment Score",
         title = paste0("Hallmark top pathways in SSc vs CTRL")) + 
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12,face = "bold"),
      axis.title = element_text(size = 14),
      legend.position = "top",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  #ggsave(filename = GSEA.H.name, plot = gsea.H.plot ,width = 12, height = 6, dpi = 300)
  print(gsea.H.plot)

path_genes <- gsea.H %>% 
                filter(pathway == "HALLMARK_INTERFERON_ALPHA_RESPONSE") %>%
                pull(leadingEdge) %>% 
                unlist()

path_genes

res2_ann <- res2_ann %>%
            mutate(
                delabel = if_else(gene %in% path_genes, gene, NA)
            )

fig.size(7,9)
fdr_threshold <- -log10(max(res2_ann$p_value[res2_ann$p_value <= 0.05], na.rm = TRUE))

volplot <- ggplot(res2_ann, aes(x = Estimate, y = -log10(adj.P.Val),label = delabel)) +
  geom_point(alpha = 0.9, 
             shape=21, 
             colour="black", 
             stroke=0.3, 
             size=2.5, 
             aes(fill=significance)) + 
  theme_minimal() +  
  labs(title = "HALLMARK_INTERFERON_ALPHA_RESPONSE genes", 
       x = "Estimate", 
       y = "-log10(Adj p-value)",
       fill=NULL,
       color=NULL) +
  geom_label_repel(data = filter(res2_ann, 
                                 !is.na(delabel)), 
                   color = "black",
                   max.overlaps = 20, 
                   force = 5, 
                   size = 3, 
                   show.legend = FALSE) +
  scale_fill_manual(values = c("Upregulated" = "#F77663", 
                               "Downregulated" = "#3D799E", 
                               "Not_significant" = "#f0f0f0")) +
  geom_hline(yintercept = fdr_threshold, linetype = "dashed", color = "blue")+
  geom_vline(xintercept = 0, linetype = 'dashed')
volplot

path_genes <- gsea.H %>% 
                filter(pathway == "HALLMARK_INTERFERON_GAMMA_RESPONSE") %>%
                pull(leadingEdge) %>% 
                unlist()

path_genes

res2_ann <- res2_ann %>%
            mutate(
                delabel = if_else(gene %in% path_genes, gene, NA)
            )

fig.size(7,9)
fdr_threshold <- -log10(max(res2_ann$p_value[res2_ann$p_value <= 0.05], na.rm = TRUE))

volplot <- ggplot(res2_ann, aes(x = Estimate, y = -log10(adj.P.Val),label = delabel)) +
  geom_point(alpha = 0.9, 
             shape=21, 
             colour="black", 
             stroke=0.3, 
             size=2.5, 
             aes(fill=significance)) + 
  theme_minimal() +  
  labs(title = "HALLMARK_INTERFERON_GAMMA_RESPONSE genes", 
       x = "Estimate", 
       y = "-log10(Adj p-value)",
       fill=NULL,
       color=NULL) +
  geom_label_repel(data = filter(res2_ann, 
                                 !is.na(delabel)), 
                   color = "black",
                   max.overlaps = 20, 
                   force = 5, 
                   size = 3, 
                   show.legend = FALSE) +
  scale_fill_manual(values = c("Upregulated" = "#F77663", 
                               "Downregulated" = "#3D799E", 
                               "Not_significant" = "#f0f0f0")) +
  geom_hline(yintercept = fdr_threshold, linetype = "dashed", color = "blue")+
  geom_vline(xintercept = 0, linetype = 'dashed')
volplot

path_genes <- gsea.H %>% 
                filter(pathway == "HALLMARK_INFLAMMATORY_RESPONSE") %>%
                pull(leadingEdge) %>% 
                unlist()

path_genes

res2_ann <- res2_ann %>%
            mutate(
                delabel = if_else(gene %in% path_genes, gene, NA)
            )

fig.size(7,9)
fdr_threshold <- -log10(max(res2_ann$p_value[res2_ann$p_value <= 0.05], na.rm = TRUE))

volplot <- ggplot(res2_ann, aes(x = Estimate, y = -log10(adj.P.Val),label = delabel)) +
  geom_point(alpha = 0.9, 
             shape=21, 
             colour="black", 
             stroke=0.3, 
             size=2.5, 
             aes(fill=significance)) + 
  theme_minimal() +  
  labs(title = "HALLMARK_INFLAMMATORY_RESPONSE genes", 
       x = "Estimate", 
       y = "-log10(Adj p-value)",
       fill=NULL,
       color=NULL) +
  geom_label_repel(data = filter(res2_ann, 
                                 !is.na(delabel)), 
                   color = "black",
                   max.overlaps = 20, 
                   force = 5, 
                   size = 3, 
                   show.legend = FALSE) +
  scale_fill_manual(values = c("Upregulated" = "#F77663", 
                               "Downregulated" = "#3D799E", 
                               "Not_significant" = "#f0f0f0")) +
  geom_hline(yintercept = fdr_threshold, linetype = "dashed", color = "blue")+
  geom_vline(xintercept = 0, linetype = 'dashed')
volplot

gsea.res <- rbind(gsea.H, gsea.c2)

fgseaResTidy <- gsea.res %>%
  as_tibble() %>%
  arrange(desc(NES))

# Show in a nice table:
fgseaResTidy %>% 
  dplyr::select( -ES) %>% 
  arrange(padj) %>% 
  DT::datatable()

top_pos <- gsea.res %>%
  filter(padj < 0.1) %>%
  filter(NES > 0) %>%
  arrange(desc(NES)) %>%
  slice_head(n = 20)

top_neg <- gsea.res %>%
  filter(padj < 0.1) %>%
  filter(NES < 0) %>%
  arrange(NES) %>%
  slice_head(n = 20)

top_combined <- bind_rows(top_pos, top_neg)

fig.size(12,14)
gsea.c2.plot <- ggplot(top_combined, aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill = ifelse(NES > 0, "NES > 0", "NES < 0")), alpha = 0.8,width = 0.75,colour="black") +
    scale_fill_manual(values = c("NES > 0" = "#FF7F50", "NES < 0" = "#6495ED"), name = "FDR 10%") +
    coord_flip() +
    labs(x = "", y = "Normalized Enrichment Score",
         title = "GSEA top pathways in SSc vs CTRL") + 
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 16),
      axis.text.y = element_text(size = 12,face = "bold"),
      axis.title = element_text(size = 16, face = "bold", hjust = 1),
      plot.title = element_text(hjust = 2, size = 16, face = "bold"),
      legend.position = "top",
      legend.title = element_text(size = 16, hjust = 4),
      legend.text = element_text(size = 16, hjust = 4),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )

gsea.c2.plot

ggsave(gsea.c2.plot,filename = "../results/plots/20260203_SSc_CD4_PseudoBulk_GSEA_combine.png",width = 14,height = 12, dpi = 300)

path_genes <- gsea.H %>% 
                filter(pathway == "HALLMARK_INTERFERON_ALPHA_RESPONSE") %>%
                pull(leadingEdge) %>% 
                unlist()

exprs_norm %>% str

#exprs_norm <- counts_store

genes_present <- path_genes[path_genes %in% rownames(exprs_norm)]
genes_missing <- path_genes[!path_genes %in% rownames(exprs_norm)]

length(genes_present)
length(genes_missing)
genes_missing

exp_heat <- exprs_norm[path_genes,]
exp_heat <- as.data.frame(t(exp_heat))

exp_heat %>% dim

exp_heat$sample <- as.character(meta_ssc$Sample)
exp_heat %>% dim
exp_heat %>% head

exp_ave <- aggregate(exp_heat[, 1:(ncol(exp_heat)-1)], list(exp_heat$sample), mean)

exp_ave %>% head

colnames(exp_ave)[1] <- "sample"
exp_ave <- as.data.frame(t(exp_ave))
colnames(exp_ave) <- as.character(t(exp_ave[1,]))
exp_ave <- exp_ave[-1,]
row_names <- rownames(exp_ave)

exp_ave <- mutate_all(exp_ave, function(x) as.numeric(as.character(x)))
rownames(exp_ave) <- row_names
exp_ave <- as.matrix(exp_ave)
scale_rows <- function(x) t(scale(t(x)))
exp_ave_scale <- scale_rows(exp_ave) # Z-score
exp_ave_scale[exp_ave_scale > 2] <- 2
exp_ave_scale[exp_ave_scale < -2] <- -2
dim(exp_ave_scale)

subset <- which(rownames(exp_ave_scale) %in% rownames(exp_ave_scale)) 
labels = rownames(exp_ave_scale)[subset]

row_annot <- data.frame(Link = labels) 
rownames(row_annot) <- rownames(exp_ave_scale)

fig.size(12,8)
pheatmap(
  exp_ave_scale,
  
  # clustering
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  
  # nombres de filas/columnas
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_col = 10,
  angle_col = 45,
  
  # estética
  treeheight_row = 0,   # equivalente a show_row_dend = FALSE
  treeheight_col = 0,   # equivalente a show_column_dend = FALSE
  
  fontsize = 12,
  border_color = NA
)

IFNa_genes <- gsea.H %>% 
                filter(pathway == "HALLMARK_INTERFERON_ALPHA_RESPONSE") %>%
                pull(leadingEdge) %>% 
                unlist()

IFNg_genes <- gsea.H %>% 
                filter(pathway == "HALLMARK_INTERFERON_GAMMA_RESPONSE") %>%
                pull(leadingEdge) %>% 
                unlist()

infla_genes <- gsea.H %>% 
                filter(pathway == "HALLMARK_INFLAMMATORY_RESPONSE") %>%
                pull(leadingEdge) %>% 
                unlist()

genes <- c(IFNa_genes,IFNg_genes,infla_genes) %>% unique

exp_heat <- exprs_norm[genes,]
exp_heat <- as.data.frame(t(exp_heat))
genes %>% length
exp_heat %>% dim
exp_heat$sample <- as.character(meta_ssc$Sample)

exp_ave <- aggregate(exp_heat[, 1:(ncol(exp_heat)-1)], list(exp_heat$sample), mean)

#exp_ave %>% head

colnames(exp_ave)[1] <- "sample"
exp_ave <- as.data.frame(t(exp_ave))
colnames(exp_ave) <- as.character(t(exp_ave[1,]))
exp_ave <- exp_ave[-1,]
row_names <- rownames(exp_ave)

exp_ave <- mutate_all(exp_ave, function(x) as.numeric(as.character(x)))

rownames(exp_ave) <- row_names
exp_ave <- as.matrix(exp_ave)
scale_rows <- function(x) t(scale(t(x)))
exp_ave_scale <- scale_rows(exp_ave) # Z-score
exp_ave_scale[exp_ave_scale > 2] <- 2
exp_ave_scale[exp_ave_scale < -2] <- -2
dim(exp_ave_scale)

ifna_score <- colMeans(exp_ave_scale[IFNa_genes, , drop = FALSE])
ifng_score <- colMeans(exp_ave_scale[IFNg_genes, , drop = FALSE])
infla_score <- colMeans(exp_ave_scale[infla_genes, , drop = FALSE])

df_score <- data.frame(
  Sample = names(ifna_score),
  ifna_score = ifna_score,
    ifng_score = ifng_score,
    infla_score = infla_score
)
#df_score

df_modules <- left_join(df_score, meta_ssc, by = "Sample")

ifna_plot <- ggplot(df_modules, aes(x=factor(condition, levels=c("CTRL","SSc")), y=ifna_score, fill=condition)) + 
  geom_boxplot(width=0.5, alpha = 0.75)+
  geom_point(alpha = 0.9, 
             shape=21, 
             colour="black", 
             stroke=0.3, 
             size=2.5) +
  labs(title = "Interferon alpha score", 
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
ifng_plot <- ggplot(df_modules, aes(x=factor(condition, levels=c("CTRL","SSc")), y=ifng_score, fill=condition)) + 
  geom_boxplot(width=0.5, alpha = 0.75)+
  geom_point(alpha = 0.9, 
             shape=21, 
             colour="black", 
             stroke=0.3, 
             size=2.5) +
  labs(title = "Interferon gamma score", 
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
infla_plot <- ggplot(df_modules, aes(x=factor(condition, levels=c("CTRL","SSc")), y=infla_score, fill=condition)) + 
  geom_boxplot(width=0.5, alpha = 0.75)+
  geom_point(alpha = 0.9, 
             shape=21, 
             colour="black", 
             stroke=0.3, 
             size=2.5) +
  labs(title = "Inflammatory response score", 
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

fig.size(4,14)
ifna_plot+ifng_plot+infla_plot

counts_con <- readRDS("2026-01-16_SSc_CD4_pseudobulk_counts_normalized.rds")

counts_con %>% str

# Remove NAs and set row names
counts_df <- as.data.frame(as.matrix(counts_con))

head(counts_df)

net <- decoupleR::get_collectri(organism = 'human', 
                                split_complexes = FALSE)

net %>% head

# Run ulm
sample_acts <- decoupleR::run_ulm(mat = counts_df, 
                                  net = net, 
                                  .source = 'source', 
                                  .target = 'target',
                                  .mor = 'mor', 
                                  minsize = 5)

sample_acts %>% head

sample_acts %>% tail

n_tfs <- 40

# Transform to wide matrix
sample_acts_mat <- sample_acts %>%
                   tidyr::pivot_wider(id_cols = 'condition', 
                                      names_from = 'source',
                                      values_from = 'score') %>%
                   tibble::column_to_rownames('condition') %>%
                   as.matrix()

# Get top tfs with more variable means across clusters
tfs <- sample_acts %>%
       dplyr::group_by(source) %>%
       dplyr::summarise(std = stats::sd(score)) %>%
       dplyr::arrange(-abs(std)) %>%
       head(n_tfs) %>%
       dplyr::pull(source)

sample_acts_mat <- sample_acts_mat[,tfs]

# Scale per sample
sample_acts_mat <- scale(sample_acts_mat)

# Choose color palette
colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
colors.use <- grDevices::colorRampPalette(colors = colors)(100)

my_breaks <- c(seq(-2, 0, length.out = ceiling(100 / 2) + 1),
               seq(0.05, 2, length.out = floor(100 / 2)))



sample_acts_mat %>% head

fig.size(8,16)
# Plot
pheatmap::pheatmap(mat = sample_acts_mat,
                   color = colors.use,
                   border_color = "white",
                   #cellwidth = 20,
                   #cellheight = 20,
                   cluster_rows = FALSE,
                   cluster_cols = TRUE,
                   treeheight_row = 20,
                   treeheight_col = 20)

res2_ann %>% dim
res2_ann %>% str

deg <- res2_ann %>%
       dplyr::select(gene, Estimate, t_value, p_value) %>% 
       dplyr::filter(!is.na(t_value)) %>% 
       tibble::column_to_rownames(var = "gene") %>%
       as.matrix()

head(deg)

# Run ulm
contrast_acts <- decoupleR::run_ulm(mat = deg[, 't_value', drop = FALSE], 
                                    net = net, 
                                    .source = 'source', 
                                    .target = 'target',
                                    .mor='mor', 
                                    minsize = 5)

contrast_acts %>% head

# Filter top TFs in both signs
f_contrast_acts <- contrast_acts %>%
                   dplyr::mutate(rnk = NA)

msk <- f_contrast_acts$score > 0

f_contrast_acts[msk, 'rnk'] <- rank(-f_contrast_acts[msk, 'score'])
f_contrast_acts[!msk, 'rnk'] <- rank(-abs(f_contrast_acts[!msk, 'score']))

tfs <- f_contrast_acts %>%
       dplyr::arrange(rnk) %>%
       head(n_tfs) %>%
       dplyr::pull(source)

f_contrast_acts <- f_contrast_acts %>%
                   filter(source %in% tfs)

f_contrast_acts %>% dim
f_contrast_acts %>% head

fig.size(6,9)
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
     xlab("TFs")

p

#ggsave(p,filename = "../results/plots/20260203_SSc_CD4_PseudoBulk_TFact.png",width = 8,height = 6, dpi = 300)

#write.csv(f_contrast_acts,"../results/GSEA_TFactivity/2026-01-31_SSc_CD4_GSEA_pseudobulk_TFact.csv")

stat1_targets <- net %>%
  filter(source == "STAT1")

stat1_targets %>% str

TF_genes <- res2_ann %>% 
            filter(res2_ann$gene %in% stat1_targets$target, res2_ann$significance!='Not_significant') %>%
            select(gene) %>%
            as.vector()

TF_genes$gene

res2_ann <- res2_ann %>%
            mutate(
                delabel = if_else(gene %in% TF_genes$gene, gene, NA)
            )

fig.size(7,9)
fdr_threshold <- -log10(max(res2_ann$p_value[res2_ann$p_value <= 0.05], na.rm = TRUE))

volplot <- ggplot(res2_ann, aes(x = Estimate, y = -log10(adj.P.Val),label = delabel)) +
  geom_point(alpha = 0.9, 
             shape=21, 
             colour="black", 
             stroke=0.3, 
             size=2.5, 
             aes(fill=significance)) + 
  scale_color_manual(values = c("Upregulated" = "#7FB77E", 
                                "Downregulated" = "#9730AF", "Not_significant" = "#bdbdbd")) +  
  theme_minimal() +  
  labs(title = "", 
       x = "Estimate", 
       y = "-log10(Adj p-value)",
       fill=NULL,
       color=NULL) +
  geom_label_repel(data = filter(res2_ann, 
                                 !is.na(delabel)), 
                   color = "black",
                   max.overlaps = 20, 
                   force = 5, 
                   size = 3, 
                   show.legend = FALSE) +
  scale_fill_manual(values = c("Upregulated" = "#7FB77E", 
                               "Downregulated" = "#da5cf7", 
                               "Not_significant" = "#f0f0f0")) +
  geom_hline(yintercept = fdr_threshold, linetype = "dashed", color = "blue")+
  geom_vline(xintercept = 0, linetype = 'dashed')
volplot

g_targets <- net %>%
  filter(source == "BMAL1")

g_targets %>% str

TF_genes <- res2_ann %>% 
            filter(res2_ann$gene %in% g_targets$target, res2_ann$significance!='Not_significant') %>%
            select(gene) %>%
            as.vector()

TF_genes$gene

res2_ann <- res2_ann %>%
            mutate(
                delabel = if_else(gene %in% TF_genes$gene, gene, NA)
            )

fig.size(7,9)
fdr_threshold <- -log10(max(res2_ann$p_value[res2_ann$p_value <= 0.05], na.rm = TRUE))

volplot <- ggplot(res2_ann, aes(x = Estimate, y = -log10(adj.P.Val),label = delabel)) +
  geom_point(alpha = 0.9, 
             shape=21, 
             colour="black", 
             stroke=0.3, 
             size=2.5, 
             aes(fill=significance)) + 
  scale_color_manual(values = c("Upregulated" = "#7FB77E", 
                                "Downregulated" = "#9730AF", "Not_significant" = "#bdbdbd")) +  
  theme_minimal() +  
  labs(title = "BMAL1 related genes in SSc vs CTRL", 
       x = "Estimate", 
       y = "-log10(Adj p-value)",
       fill=NULL,
       color=NULL) +
  geom_label_repel(data = filter(res2_ann, 
                                 !is.na(delabel)), 
                   color = "black",
                   max.overlaps = 20, 
                   force = 5, 
                   size = 3, 
                   show.legend = FALSE) +
  scale_fill_manual(values = c("Upregulated" = "#7FB77E", 
                               "Downregulated" = "#da5cf7", 
                               "Not_significant" = "#f0f0f0")) +
  geom_hline(yintercept = fdr_threshold, linetype = "dashed", color = "blue")+
  geom_vline(xintercept = 0, linetype = 'dashed')
volplot

# Ranking the genes for GSEA
ranking <- setNames(contrast_acts$score, contrast_acts$source) # asign to each gene its OR

ranking %>% str

c2_gene_sets = msigdbr(species = "human", category = 'C2')

react <- c2_gene_sets %>%
  filter(gs_subcollection == "CP:REACTOME") %>%
  select(gs_name, gene_symbol) %>% 
  group_by(gs_name) %>% 
  summarise(all.genes = list(unique(gene_symbol))) %>% 
  deframe()

C2.ensembl.ls <- c2_gene_sets %>% 
  select(gs_name, gene_symbol) %>% 
  group_by(gs_name) %>% 
  summarise(all.genes = list(unique(gene_symbol))) %>% 
  deframe()

gsea.c2 <- fgseaMultilevel(pathways= react,
                           stats=ranking,  
                           minSize = 5, 
                           maxSize = 500)
head(gsea.c2[order(pval), ])

fgseaResTidy <- gsea.c2 %>%
  as_tibble() %>%
  arrange(desc(NES))

# Show in a nice table:
fgseaResTidy %>% 
  dplyr::select( -ES) %>% 
  arrange(padj) %>% 
  DT::datatable()

top_pos <- gsea.c2 %>%
  filter(padj < 0.1) %>%
  filter(NES > 0) %>%
  arrange(desc(NES)) %>%
  slice_head(n = 10)

top_neg <- gsea.c2 %>%
  filter(padj < 0.1) %>%
  filter(NES < 0) %>%
  arrange(NES) %>%
  slice_head(n = 10)

top_combined <- bind_rows(top_pos, top_neg)

fig.size(8,13)
gsea.c2.plot <- ggplot(top_combined, aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill = ifelse(NES > 0, "NES > 0", "NES < 0")), alpha = 0.8,width = 0.75,colour="black") +
    scale_fill_manual(values = c("NES > 0" = "#FF7F50", "NES < 0" = "#6495ED"), name = "FDR 10%") +
    coord_flip() +
    labs(x = "Pathway", y = "Normalized Enrichment Score",
         title = paste0("REACTOME top pathways in SSc vs CTRL")) + 
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12,face = "bold"),
      axis.title = element_text(size = 14),
      legend.position = "top",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )

gsea.c2.plot

fig.size(6,10)
plotEnrichment(react[["REACTOME_INTERFERON_ALPHA_BETA_SIGNALING"]], ranking) +
  labs(
    title = "IFNa and IFNb signaling"
  )

fig.size(6,10)
plotEnrichment(react[["REACTOME_PYROPTOSIS"]], ranking) +
  labs(
    title = "Pyroptosis"
  )

fig.size(6,10)
plotEnrichment(react[["REACTOME_NGF_STIMULATED_TRANSCRIPTION"]], ranking) +
  labs(
    title = ""
  )

gsea.H <- fgseaMultilevel(pathways = H.ensembl.ls, 
                           stats=ranking,  
                           minSize = 10, 
                           maxSize = 500)


fgseaResTidy <- gsea.H %>%
  as_tibble() %>%
  arrange(desc(NES))

# Show in a nice table:
fgseaResTidy %>% 
  dplyr::select(-ES) %>% 
  arrange(padj) %>% 
  DT::datatable()

fig.size(6,10)
plotEnrichment(H.ensembl.ls[["HALLMARK_INFLAMMATORY_RESPONSE"]], ranking) +
  labs(
    title = "Inflammatory response"
  )

plotEnrichment(H.ensembl.ls[["HALLMARK_TNFA_SIGNALING_VIA_NFKB"]], ranking) +
  labs(
    title = "TNFa via NKFB"
  )

plotEnrichment(H.ensembl.ls[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]], ranking) +
  labs(
    title = "IFNg"
  )

top_pos.H <- gsea.H %>%
  filter(padj < 0.1) %>%
  filter(NES > 0) %>%
  arrange(desc(NES)) %>%
  slice_head(n = 10)

top_neg.H <- gsea.H %>%
  filter(padj < 0.1) %>%
  filter(NES < 0) %>%
  arrange(NES) %>%
  slice_head(n = 10)

top_combined.H <- bind_rows(top_pos.H, top_neg.H)

fig.size(8,10)
gsea.H.plot <- ggplot(top_combined.H, aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill = ifelse(NES > 0, "NES > 0", "NES < 0")), alpha = 0.8,width = 0.75,colour="black") +
    scale_fill_manual(values = c("NES > 0" = "#FF7F50", "NES < 0" = "#6495ED"), name = "FDR 10%") +
    coord_flip() +
    labs(x = "Pathway", y = "Normalized Enrichment Score",
         title = paste0("Hallmark top pathways in SSc vs CTRL")) + 
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12,face = "bold"),
      axis.title = element_text(size = 14),
      legend.position = "top",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  #ggsave(filename = GSEA.H.name, plot = gsea.H.plot ,width = 12, height = 6, dpi = 300)
  print(gsea.H.plot)

gene_score <- function(genes, meta, exprs_norm, score_name){

  ## -------------------------------
  ## 1. Filtrar genes presentes
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
  exp_heat <- as.data.frame(t(exp_heat))

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

act_genes <- c(
  "CD69",
  "IL2RA",
  "ICOS",
  "TNFRSF4",
  "TNFRSF9",
  "FOS",
  "JUNB"
)

genes <- act_genes

genes_present <- genes[genes %in% rownames(exprs_norm)]
genes_missing <- genes[!genes %in% rownames(exprs_norm)]

length(genes_present)
length(genes_missing)
genes_missing

exp_heat <- exprs_norm[genes,]
exp_heat <- as.data.frame(t(exp_heat))
genes %>% length
exp_heat %>% dim
exp_heat$sample <- as.character(meta_ssc$Sample)

exp_ave <- aggregate(exp_heat[, 1:(ncol(exp_heat)-1)], list(exp_heat$sample), mean)

colnames(exp_ave)[1] <- "sample"
exp_ave <- as.data.frame(t(exp_ave))
colnames(exp_ave) <- as.character(t(exp_ave[1,]))
exp_ave <- exp_ave[-1,]
row_names <- rownames(exp_ave)

exp_ave <- mutate_all(exp_ave, function(x) as.numeric(as.character(x)))

rownames(exp_ave) <- row_names
exp_ave <- as.matrix(exp_ave)
scale_rows <- function(x) t(scale(t(x)))
exp_ave_scale <- scale_rows(exp_ave) # Z-score
exp_ave_scale[exp_ave_scale > 2] <- 2
exp_ave_scale[exp_ave_scale < -2] <- -2
dim(exp_ave_scale)

score <- colMeans(exp_ave_scale[genes, , drop = FALSE])

df_score <- data.frame(
  Sample = names(score),
  score = score
)
#df_score

df_modules <- left_join(df_score, meta_ssc, by = "Sample")

fig.size(6,6)
act_modulePlot <- ggplot(df_modules, aes(x=factor(condition, levels=c("CTRL","SSc")), y=score, fill=condition)) + 
  geom_boxplot(width=0.5, alpha = 0.75)+
  geom_point(alpha = 0.9, 
             shape=21, 
             colour="black", 
             stroke=0.3, 
             size=2.5) +
  labs(title = "Activation score", 
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
act_modulePlot

res2_ann %>% colnames

res2_exh <- res2_ann

res2_exh$exh <- ifelse(res2_exh$gene %in% genes, "Activation", "Other")
res2_exh$delabel <- NA
res2_exh$delabel[res2_exh$exh != "Other"] <- res2_exh$gene[res2_exh$exh != "Other"]

fig.size(8,10)
volplot_act <- ggplot(res2_exh, aes(x = Estimate, y = -log10(p_value))) +

  geom_point(
    data = subset(res2_exh, exh == "Other"),
    shape = 21,
    fill = "#f0f0f0",
    colour = "black",
    stroke = 0.2,
    size = 2,
    alpha = 0.6
  ) +

  geom_point(
    data = subset(res2_exh, exh == "Activation"),
    aes(fill = exh),
    shape = 21,
    colour = "black",
    stroke = 0.4,
    size = 3.5
  ) +

  theme_minimal() +
  
  theme(legend.position = "") +

  geom_label_repel(
    data = subset(res2_exh, gene %in% genes),
    aes(label = gene),
    color = "black",
    size = 3.5,
    max.overlaps = Inf,
    force = 5,
    show.legend = FALSE
  ) +

  scale_fill_manual(values = c("Activation" = "#31a354")) +

  geom_hline(yintercept = fdr_threshold, linetype = "dashed", color = "blue")+
  geom_vline(xintercept = 0, linetype = "dashed") +

  labs(
    x = "Estimate",
    y = "-log10(p-value)",
    fill = ""
  )

fig.size(4,6)
volplot_act

fig.size(4,10)
act_modulePlot+volplot_act

exha <- c("PDCD1","HAVCR2","LAG3","TOX","TOX2","NR4A1","NR4A2","NR4A3","PRDM1","ENTPD1")

genes <- exha

#exprs_norm <- counts_store

genes_present <- genes[genes %in% rownames(exprs_norm)]
genes_missing <- genes[!genes %in% rownames(exprs_norm)]

length(genes_present)
length(genes_missing)
genes_missing

exp_heat <- exprs_norm[genes,]
exp_heat <- as.data.frame(t(exp_heat))
genes %>% length
exp_heat %>% dim
exp_heat$sample <- as.character(meta_ssc$Sample)

exp_ave <- aggregate(exp_heat[, 1:(ncol(exp_heat)-1)], list(exp_heat$sample), mean)

colnames(exp_ave)[1] <- "sample"
exp_ave <- as.data.frame(t(exp_ave))
colnames(exp_ave) <- as.character(t(exp_ave[1,]))
exp_ave <- exp_ave[-1,]
row_names <- rownames(exp_ave)

exp_ave <- mutate_all(exp_ave, function(x) as.numeric(as.character(x)))

rownames(exp_ave) <- row_names
exp_ave <- as.matrix(exp_ave)
scale_rows <- function(x) t(scale(t(x)))
exp_ave_scale <- scale_rows(exp_ave) # Z-score
exp_ave_scale[exp_ave_scale > 2] <- 2
exp_ave_scale[exp_ave_scale < -2] <- -2
dim(exp_ave_scale)

score <- colMeans(exp_ave_scale[genes, , drop = FALSE])

df_score <- data.frame(
  Sample = names(score),
  score = score
)
#df_score

df_modules <- left_join(df_score, meta_ssc, by = "Sample")

fig.size(6,6)
ex_modulePlot <- ggplot(df_modules, aes(x=factor(condition, levels=c("CTRL","SSc")), y=score, fill=condition)) + 
  geom_boxplot(width=0.5, alpha = 0.75)+
  geom_point(alpha = 0.9, 
             shape=21, 
             colour="black", 
             stroke=0.3, 
             size=2.5) +
  labs(title = "Exhaustion score", 
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
ex_modulePlot

res2_ann %>% colnames

res2_exh <- res2_ann

res2_exh$exh <- ifelse(res2_exh$gene %in% genes, "Exhaustion", "Other")
res2_exh$delabel <- NA
res2_exh$delabel[res2_exh$exh != "Other"] <- res2_exh$gene[res2_exh$exh != "Other"]

fig.size(8,10)
volplot_exh <- ggplot(res2_exh, aes(x = Estimate, y = -log10(p_value))) +

  geom_point(
    data = subset(res2_exh, exh == "Other"),
    shape = 21,
    fill = "#f0f0f0",
    colour = "black",
    stroke = 0.2,
    size = 2,
    alpha = 0.6
  ) +

  geom_point(
    data = subset(res2_exh, exh == "Exhaustion"),
    aes(fill = exh),
    shape = 21,
    colour = "black",
    stroke = 0.4,
    size = 3.5
  ) +

  theme_minimal() +
  
  theme(legend.position = "") +

  geom_label_repel(
    data = subset(res2_exh, gene %in% genes),
    aes(label = gene),
    color = "black",
    size = 3.5,
    max.overlaps = Inf,
    force = 5,
    show.legend = FALSE
  ) +

  scale_fill_manual(values = c("Exhaustion" = "#31a354")) +

  geom_hline(yintercept = fdr_threshold, linetype = "dashed", color = "blue")+
  geom_vline(xintercept = 0, linetype = "dashed") +

  labs(
    x = "Estimate",
    y = "-log10(p-value)",
    fill = ""
  )

fig.size(4,6)
volplot_exh

fig.size(4,10)
ex_modulePlot+volplot_exh

IFNa_genes <- gsea.H %>% 
                filter(pathway == "HALLMARK_INTERFERON_ALPHA_RESPONSE") %>%
                pull(leadingEdge) %>% 
                unlist()

IFNg_genes <- gsea.H %>% 
                filter(pathway == "HALLMARK_INTERFERON_GAMMA_RESPONSE") %>%
                pull(leadingEdge) %>% 
                unlist()

genes <- c(IFNa_genes,IFNg_genes) %>% unique

genes

#genes <- c("IFITM1","MX1","MX2","IFITM2","IFI35","STAT1","STAT2","STAT3","ISG15","IFNG")

#exprs_norm <- counts_store

genes_present <- genes[genes %in% rownames(exprs_norm)]
genes_missing <- genes[!genes %in% rownames(exprs_norm)]

length(genes_present)
length(genes_missing)
genes_missing

exp_heat <- exprs_norm[genes,]
exp_heat <- as.data.frame(t(exp_heat))
genes %>% length
exp_heat %>% dim
exp_heat$sample <- as.character(meta_ssc$Sample)

exp_ave <- aggregate(exp_heat[, 1:(ncol(exp_heat)-1)], list(exp_heat$sample), mean)

colnames(exp_ave)[1] <- "sample"
exp_ave <- as.data.frame(t(exp_ave))
colnames(exp_ave) <- as.character(t(exp_ave[1,]))
exp_ave <- exp_ave[-1,]
row_names <- rownames(exp_ave)

exp_ave <- mutate_all(exp_ave, function(x) as.numeric(as.character(x)))

rownames(exp_ave) <- row_names
exp_ave <- as.matrix(exp_ave)
scale_rows <- function(x) t(scale(t(x)))
exp_ave_scale <- scale_rows(exp_ave) # Z-score
exp_ave_scale[exp_ave_scale > 2] <- 2
exp_ave_scale[exp_ave_scale < -2] <- -2
dim(exp_ave_scale)

score <- colMeans(exp_ave_scale[genes, , drop = FALSE])

df_score <- data.frame(
  Sample = names(score),
  score = score
)
#df_score

df_modules <- left_join(df_score, meta_ssc, by = "Sample")

fig.size(6,6)
ifn_modulePlot <- ggplot(df_modules, aes(x=factor(condition, levels=c("CTRL","SSc")), y=score, fill=condition)) + 
  geom_boxplot(width=0.5, alpha = 0.75)+
  geom_point(alpha = 0.9, 
             shape=21, 
             colour="black", 
             stroke=0.3, 
             size=2.5) +
  labs(title = "IFN score", 
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
ifn_modulePlot

res2_ann %>% colnames

res2_ifn <- res2_ann

res2_ifn$ifn <- ifelse(res2_ifn$gene %in% genes, "IFN", "Other")
res2_ifn$delabel <- NA
res2_ifn$delabel[res2_ifn$ifn != "Other"] <- res2_ifn$gene[res2_ifn$ifn != "Other"]

fig.size(8,10)
fdr_threshold <- -log10(max(res2_ann$p_value[res2_ann$adj.P.Val <= 0.05], na.rm = TRUE))
volplot <- ggplot(res2_ifn, aes(x = Estimate, y = -log10(p_value))) +

  geom_point(
    data = subset(res2_ifn, ifn == "Other"),
    shape = 21,
    fill = "#f0f0f0",
    colour = "black",
    stroke = 0.2,
    size = 2,
    alpha = 0.6
  ) +

  geom_point(
    data = subset(res2_ifn, ifn == "IFN"),
    aes(fill = ifn),
    shape = 21,
    colour = "black",
    stroke = 0.4,
    size = 3.5
  ) +

  geom_label_repel(
    data = subset(res2_ifn, gene %in% genes),
    aes(label = gene),
    color = "black",
    size = 3.5,
    max.overlaps = 20,
    force = 5,
    show.legend = FALSE
  ) +

  scale_fill_manual(values = c("IFN" = "#c51b8a")) +

  geom_hline(yintercept = fdr_threshold, linetype = "dashed", color = "blue")+
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal() +
  theme(legend.position = "") +
  labs(
    x = "Estimate",
    y = "-log10(p-value)",
    fill = ""
  )

fig.size(4,6)
volplot

fig.size(4,10)
ifn_modulePlot+volplot

HLA_cii <- c("HLA-DRA","HLA-DRB1","HLA-DQA1","HLA-DQB1","HLA-DPA1","HLA-DPB1",
             "CD69",
             "IL2RA"
            )

genes <- HLA_cii

hla_score <- gene_score(HLA_cii, meta_ssc, exprs_norm, "MHC-II_score")

df_score <- data.frame(
  Sample = names(hla_score),
  hla_score = hla_score
)

hla_plot <- moduleBoxPlot(df_modules,df_modules$hla_score,"Activation-IFN")

res2_hla <- res2_ann

res2_hla$act <- case_when(
    res2_hla$gene %in% HLA_cii ~ "Activation-IFN",
    TRUE ~ "Other"
)

#res2_hla$hla <- ifelse(res2_hla$gene %in% genes, "Activation", "Other")
res2_hla$delabel <- NA
res2_hla$delabel[res2_hla$act != "Other"] <- res2_hla$gene[res2_hla$act != "Other"]

fig.size(8,10)
MHCplot <- ggplot(res2_hla, aes(x = Estimate, y = -log10(p_value))) +

  geom_point(
    data = subset(res2_hla, act == "Other"),
    shape = 21,
    fill = "#f0f0f0",
    colour = "black",
    stroke = 0.2,
    size = 2,
    alpha = 0.6
  ) +

  geom_point(
    data = subset(res2_hla, act == "Activation-IFN"),
    aes(fill = act),
    shape = 21,
    colour = "black",
    stroke = 0.4,
    size = 3.5
  ) +

  geom_label_repel(
    data = subset(res2_hla, gene %in% genes),
    aes(label = gene),
    color = "black",
    size = 3.5,
    max.overlaps = Inf,
    force = 5,
    show.legend = FALSE
  ) +

  scale_fill_manual(values = c("Activation-IFN" = "#3182bd")) +

  geom_hline(yintercept = fdr_threshold, linetype = "dashed", color = "blue")+
  geom_vline(xintercept = 0, linetype = "dashed") +
  

  theme_minimal() +
  theme(legend.position = "") +
  
  labs(
    x = "Estimate",
    y = "-log10(p-value)",
    fill = ""
  )

fig.size(4,6)
MHCplot

fig.size(4,10)
hla_plot+MHCplot

# Modules to plot: Exhaustion, MHC and IFN
# ex_modulePlot+volplot_exh
# hla_plot+MHCplot
# ifn_modulePlot+volplot

# Combine plots plots
box_score <- wrap_plots(ex_modulePlot,hla_plot,ifn_modulePlot)
vol_score <- wrap_plots(volplot_exh,MHCplot,volplot)

fig.size(8,14)
modules_final <- wrap_plots(box_score, vol_score, ncol = 1)
modules_final

#ggsave(modules_final,filename = "../results/plots/20260204_SSc_CD4_PseudoBulk_modules_v2.png",width = 14,height = 8, dpi = 300)



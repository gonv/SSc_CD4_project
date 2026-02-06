#!/usr/bin/env Rscript

<<<<<<< HEAD
data_lib <- c("devtools", "tidyverse", "ggpubr", "ggrepel", "enrichR","matrixStats",
              "RColorBrewer", "reshape2","shiny","DT","ggpubr","DelayedArray",
              "bslib","shinyjs","ggsignif","ggdendro","cowplot","pheatmap","rstatix", "Seurat")
suppressPackageStartupMessages({
  invisible(lapply(data_lib, require, character.only = TRUE))
})
=======
# volcano plot
volcano_DEA <- function(df, tg, SCR,tgStim,tgCon,genesToRep = NULL) {
  # Initialize additional columns
	
	res2_ann <- df
  res2_ann$logP <- -log10(res2_ann$p_value)
  res2_ann$delabel <- NA

  # Generate title for the plot
  titleN <- paste(tg, "KO vs", SCR, tgStim, tgCon, sep = " ")

  # Calculate the FDR threshold for visualization
  fdr_threshold <- -log10(max(res2_ann$p_value[res2_ann$q_value <= 0.05], na.rm = TRUE))
  
  # Create the plot 
  volcano.p <- ggplot(res2_ann, aes(x = Beta, y = logP, color = Result)) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = c("Upregulated" = "#7FB77E", "Downregulated" = "#9730AF", "Not_significant" = "#bdbdbd")) +
    geom_hline(yintercept = fdr_threshold, linetype = "dashed", color = "blue") +
    ggtitle(titleN)  # Adding the title to the plot
  
  if(!is.null(genesToRep)){
    res2_ann <- res2_ann %>%
      mutate(
        delabel = ifelse(geneName %in% genesToRep, geneName, delabel)  # Label only genes in gene_vector
      )
    
    volcano.p <- volcano.p + 
      ggrepel::geom_label_repel(data = dplyr::filter(res2_ann, !is.na(delabel)),
                                aes(label = delabel),
                                ,max.overlaps = Inf,
                                color = "black")
  }
  
  return(volcano.p)
}

# heatmap

heatmap_DEA_2 <- function(df,tg, SCR, tgStim, tgCon,logCPM,genes) {
  
  gene.vec <-genes
  geneID.vec <- unlist(sapply(gene.vec, gxtotx_conver))
  geneID.vec.df <- as.data.frame(geneID.vec)
  geneID.vec.df$geneNames <- rownames(geneID.vec.df)
  rownames(geneID.vec.df) <- geneID.vec.df$geneID.vec
  
  log_filtered.list <- logCPM[as.vector(geneID.vec.df$geneID.vec),]
  log_filtered.list <- na.omit(log_filtered.list)
  
  log_filtered.list$geneID <- rownames(log_filtered.list)
  
  log_filtered.list.2 <- merge(log_filtered.list, geneID.vec.df, by.x = "geneID", by.y = "geneID.vec")
  rownames(log_filtered.list.2) <- log_filtered.list.2$geneNames
  log_filtered.list.2 <- log_filtered.list.2[,-c(1,10)]
  
  log_filtered.list.mx <- as.matrix(log_filtered.list.2)
  
  print(str(log_filtered.list.mx))
  
  heatmap_plot <- pheatmap(log_filtered.list.mx, 
                           scale = "row",  # Normaliza los datos por fila
                           clustering_distance_rows = "euclidean",
                           clustering_distance_cols = "euclidean",
                           clustering_method = "average",  # Método de clustering, puede ser "average", "single", etc.
                           color = colorRampPalette(c("blue", "white", "red"))(100),  # Gradiente de color
                           show_rownames = TRUE,
                           show_colnames = TRUE,
                           fontsize_row = 10,
                           fontsize_col = 10,
                           angle_col = 90,  # Ángulo de los nombres de las columnas
                           main = "")  # Título del heatmap
  
  return(heatmap_plot)
}


# Boxplot

boxplot_DEA <- function(tgGene,tg,SCR,tgStim,tgCon,metaSample,txi_gx){
  
  res <- tpm_filter(
    metaSample = metaSample, 
    txi_tpm = txi_gx,
    tg_name = tg, 
    SCR = "SCR", 
    stim = tgStim,
    condition = tgCon
  )
  # Asigning results
  txi_counts_filtered_byCon <- res$txi_tpm_filtered_byCon
  metaSubset <- res$metaSubset
  
  logCPM <- log2(txi_counts_filtered_byCon+1)
  
  KO_n <- gxtotx_conver(tgGene)
  
  logCPM_fil <- logCPM %>% filter(rownames(logCPM)==KO_n)
  
  df <- as.data.frame(t(logCPM_fil))
  colnames(df) <- c("KO")
  
  df$samples <- rownames(df)
  
  df_fnl <- merge(df,metaSubset,by="samples")
  
  plot.title <- paste(tg,"KO vs",SCR, sep = " ") 
  
  #ggplot(df_fnl, aes(x = geneKO, y = KO, color = geneKO)) +
   # geom_dotplot(binaxis = "y", stackdir = "center",aes(fill=Donor))+
  #  geom_boxplot(width=0.3,alpha = 0.5)+
   # ylab("log2CPM")+
    #scale_x_discrete(labels=c(tg = "KO", SCR = "SCR"))+
    #ggtitle(plot.title)
  
  tmp <- colnames(txi_counts_filtered_byCon)
  metaSubset_sorted <- metaSubset[match(tmp, metaSubset$samples), ]
  
  tmp_log_t <- as.data.frame(t(logCPM_fil))
  tmp_log_t$samples <- rownames(tmp_log_t)
  tmp_t_df <- merge(tmp_log_t,metaSubset_sorted, by="samples")
  tmp_t_df2 <- tmp_t_df %>%
    select(KO_n,geneKO,Donor)
  
  colnames(tmp_t_df2) <- c("log2CPM","geneKO","Donor")
  
  stat.test <- tmp_t_df2  %>% 
    t_test(log2CPM ~ geneKO, paired = TRUE) %>%
    add_significance()
  stat.test <- stat.test %>% add_xy_position(x = "geneKO")
  
  pairbox <- ggplot(tmp_t_df2, aes(x = geneKO, y = log2CPM, color = geneKO)) +
    geom_dotplot(binaxis = "y", stackdir = "center",aes(fill=Donor))+
    geom_boxplot(width = 0.5, alpha = 0.7) +
    geom_line(aes(group = Donor), color = "gray70", size = 0.4) +  
    scale_x_discrete(limits = c(SCR, tg)) + # Ordenar SCR e IRF5
    ylab("log2CPM") +
    xlab(tgGene) +
    theme_minimal() +
    theme(legend.position = "right") +
    stat_pvalue_manual(stat.test, label = "p.signif", tip.length = 0,label.size = 7) +
    labs(subtitle = get_test_label(stat.test, detailed = TRUE))
  
  return(pairbox)
  
}

# GSEA

GSEA <- function(df,FDR_out = 1, GSEA_analysis){
  
  # Define required packages
  data_lib <- c("devtools", "tidyverse", "ggpubr", "ggrepel", "enrichR",
                "RColorBrewer", "fgsea", "msigdbr","reshape2")
  
  # Check and install missing packages
  missing_libs <- data_lib[!data_lib %in% installed.packages()[, "Package"]]
  if (length(missing_libs) > 0) {
    install.packages(missing_libs, dependencies = TRUE)
  }
  
  # Load all packages silently
  suppressPackageStartupMessages({
    invisible(lapply(data_lib, require, character.only = TRUE))
  })
  
  # Set seed for reproducibility
  set.seed(42)
  
  gene_list.filtered <- df %>%
    filter(FDR < FDR_out) # we can filter genes by q_value or FDR and change the threshold here
  gene_list_raw <- df
  gene_list <- gene_list.filtered
  
  # Ranking genes for the analysis
  ranking <- setNames(gene_list$Beta, gene_list$geneName)
  
  if(GSEA_analysis=="C2"){
    
    #download the category db
    c2_gene_sets = msigdbr(species = "human", category = 'C2')
 
    C2.ensembl.ls <- c2_gene_sets %>% 
      select(gs_name, gene_symbol) %>% 
      group_by(gs_name) %>% 
      summarise(all.genes = list(unique(gene_symbol))) %>% 
      deframe()
    # performing GSEA analysis
    gsea.c2 <- fgseaSimple(pathways = C2.ensembl.ls,
                           stats = ranking,
                           nperm=10000)
    fgseaResTidy <- gsea.c2 %>%
      as_tibble() %>%
      arrange(desc(NES))
    
    # Show in a nice table:
    fgseaResTidy %>% 
      dplyr::select(-leadingEdge, -ES) %>% 
      arrange(padj) %>% 
      DT::datatable()
    
    # boton to show the table
    
    # plotting section
    
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
    
    if(!is_empty(top_combined$pval)){
      GSEA.C2.name <- paste("~/Documents/SplicingFolder/plots/",tg,"KOvs",SCR,"_",tgCon,tgStim,"_GSEA_C2_topcomb.png",sep = "")
      ftitle.c2 <- paste("Curate gene sets (C2) enrichment analysis of",tg,"vs",SCR,"in",tgStim,"and",tgCon,sep=" ")
      gsea.c2.plot <- ggplot(top_combined, aes(reorder(pathway, NES), NES)) +
        geom_col(aes(fill = ifelse(NES > 0, "NES > 0", "NES < 0")), alpha = 0.8) +
        scale_fill_manual(values = c("NES > 0" = "#FF7F50", "NES < 0" = "#6495ED"), name = "FDR 5%") +
        coord_flip() +
        labs(x = "Pathway", y = "Normalized Enrichment Score",
             title = ftitle.c2) + 
        theme_bw() +
        theme(
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 10),
          axis.title = element_text(size = 14),
          legend.position = "top",
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10),
          plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
        )
      ggsave(filename = GSEA.C2.name, plot = gsea.c2.plot ,width = 12, height = 6, dpi = 300)
      print(gsea.c2.plot)
    }else{
      print("empty table")
    }
    
    
    
  }else if(GSEA_analysis=="C5"){
    
  }else if(GSEA_analysis=="Hallmark"){
    
  }
  
  return()
}

GSEA_table <- function(){

}
>>>>>>> f40c6eb162999aea48b4d3f321b5e2bce1953d15

#Feature plot

FeaPlot <- function(geneName, meta, expres_norm) {
  gene_to_plot <- geneName  
  gene_expr <- expres_norm[gene_to_plot, ]
  
  # Recorta valores extremos para mejorar contraste
  max_cap <- quantile(gene_expr[gene_expr > 0], 0.99, na.rm = TRUE)
  gene_expr[gene_expr > max_cap] <- max_cap
  
  # Reescala entre 0 y 1 para visualización
  gene_expr_scaled <- scales::rescale(gene_expr, to = c(0, 1))
  
  meta$gene_expr <- gene_expr_scaled[meta$X]
  
  # Ordena el metadata para que los puntos con menor expresión (gris)
  # se dibujen primero, y los de mayor expresión (rojo) al final
  meta <- meta[order(meta$gene_expr, decreasing = FALSE), ]
  
  # Construye el plot
  feaPlot <- ggplot(meta, aes(x = UMAP_1, y = UMAP_2, color = gene_expr)) +
    geom_point(size = 2) +
    scale_color_gradientn(colors = c("lightgrey", "red")) +
    labs(title = paste("Expression of", gene_to_plot),
         color = "Scale expression") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank()
    )+
    coord_fixed()
  
  return(feaPlot)
}
<<<<<<< HEAD

# Function 1: Individual-level boxplot
plot_boxplot_expr_individual <- function(gene, group_var, exprs_norm, meta) {
  meta$clusters <- as.character(meta$R09)
  gene_expr <- as.numeric(Matrix::as.matrix(exprs_norm[gene, , drop = FALSE]))
  names(gene_expr) <- colnames(exprs_norm)
  
  df <- meta %>%
    mutate(expression = gene_expr[colnames(exprs_norm)]) %>%
    group_by(Sample, .data[[group_var]]) %>%
    summarise(mean_expr = mean(expression, na.rm = TRUE), .groups = "drop")
  
  custom_colors <- c(
    "CTRL" = "#D8BFD8", "SSc" = "#C3E5C3",
    "Tn"= "#a9dfbb", "Tcm"= "#6cd0c4", "Treg"= "#fc6daf",
    "Tem"= "#3faccd", "CTL"= "#ff0000","Temra"= "#ff0000",
    "0" = "#CDB4DB", "1" = "#F78181", "2" = "#82CFFD", "3" = "#8FCFD1",
    "4" = "#AEB6BF", "5" = "#F7ACB4", "6" = "#FFB347", "7" = "#7FFFD4",
    "8" = "#F4D03F", "9" = "#CDB4DB", "10" = "#FFA07A", "11" = "#2ECC71"
  )
  
  if (group_var == "cell_type") {
    df[[group_var]] <- factor(df[[group_var]], levels = c("Tn", "Tcm", "Tem", "CTL", "Treg"))
  } else if (group_var == "clusters") {
    df[[group_var]] <- factor(df[[group_var]], levels = as.character(0:11))
  } else if (group_var == "condition") {
    df[[group_var]] <- factor(df[[group_var]], levels = c("CTRL", "SSc"))
  }
  
  ggplot(df, aes(x = .data[[group_var]],
                 y = mean_expr,
                 fill = .data[[group_var]])) +
    geom_boxplot(width = 0.5, alpha = 0.75, outlier.shape = NA) +
    geom_point(alpha = 0.9,
               shape = 21,
               colour = "black",
               stroke = 0.3,
               size = 2.5,
               position = position_jitter(width = 0.15)) +
    labs(title = paste("Expression of", gene),
         x = "",
         y = "Mean expression per sample",
         fill = "") +
    scale_fill_manual(values = custom_colors, na.translate = FALSE) +
    theme_bw() +
    theme(legend.position = 'none',
          axis.text.x = element_text(size = 12, face = "bold"),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 20),
          title = element_text(size = 20),
          plot.margin = unit(c(1, 1, 1, 1), "mm"))
}

# Function 2: Grouped by variable (e.g., R09 or cell_type)
plot_boxplot_expr_grouped <- function(gene, group_var, exprs_norm, meta) {
  meta$clusters <- as.character(meta$R09)
  gene_expr <- as.numeric(Matrix::as.matrix(exprs_norm[gene, , drop = FALSE]))
  names(gene_expr) <- colnames(exprs_norm)
  
  df <- meta %>%
    mutate(expression = gene_expr[colnames(exprs_norm)]) %>%
    select(Sample, condition, clusters, cell_type, expression)
  
  group_vars <- c("Sample", group_var)
  if (group_var != "condition") group_vars <- c(group_vars, "condition")
  
  df_summary <- df %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(mean_expr = mean(expression, na.rm = TRUE), .groups = "drop")
  
  df_summary$facet_group <- factor(df_summary[[group_var]],
                                   levels = if (group_var == "cell_type") {
                                     c("Tn", "Tcm", "Tem", "CTL", "Treg")
                                   } else {
                                     as.character(0:11)
                                   })
  
  custom_colors <- c(
    "CTRL" = "#D8BFD8", "SSc" = "#C3E5C3",
    "Tn"= "#a9dfbb", "Tcm"= "#6cd0c4", "Treg"= "#fc6daf",
    "Tem"= "#3faccd", "CTL"= "#ff0000",
    "0" = "#CDB4DB", "1" = "#F78181", "2" = "#82CFFD", "3" = "#8FCFD1",
    "4" = "#AEB6BF", "5" = "#F7ACB4", "6" = "#FFB347", "7" = "#7FFFD4",
    "8" = "#F4D03F", "9" = "#CDB4DB", "10" = "#FFA07A", "11" = "#2ECC71"
  )
  
  p <- ggplot(df_summary, aes(x = condition, y = mean_expr, fill = .data[[group_var]])) +
    geom_boxplot(width = 0.5, alpha = 0.75, outlier.shape = NA) +
    geom_point(alpha = 0.9, shape = 21, colour = "black", stroke = 0.3,
               size = 2.5, position = position_jitter(width = 0.15)) +
    labs(title = paste("Expression of", gene, "(grouped by", group_var, ")"),
         x = group_var, y = "Mean expression per sample", fill = group_var) +
    scale_fill_manual(values = custom_colors, na.translate = FALSE) +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 12, face = "bold", angle = 45, hjust = 1),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 20),
          title = element_text(size = 18),
          plot.margin = unit(c(1, 1, 1, 1), "mm"))
  
  if (group_var != "condition") {
    p <- p + facet_wrap(~facet_group)
  }
  
  p
}
=======
>>>>>>> f40c6eb162999aea48b4d3f321b5e2bce1953d15

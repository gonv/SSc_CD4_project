#!/usr/bin/env Rscript

data_lib <- c("devtools", "tidyverse", "ggpubr", "ggrepel", "enrichR","matrixStats",
              "RColorBrewer", "reshape2","shiny","DT","ggpubr","DelayedArray",
              "bslib","shinyjs","ggsignif","ggdendro","cowplot","pheatmap","rstatix", "Seurat")
suppressPackageStartupMessages({
  invisible(lapply(data_lib, require, character.only = TRUE))
})

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
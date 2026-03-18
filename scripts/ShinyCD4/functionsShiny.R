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

# Function to plot gene expression on UMAP with optional separation by condition
plot_umap_expression <- function(geneName, meta, expres_norm, split_mode = "Combined") {
  
  gene_to_plot <- geneName
  
  # Extraer expresión del gen
  if (!gene_to_plot %in% rownames(expres_norm)) {
    warning("Gene not found in expression matrix.")
    return(ggplot() + ggtitle("Gene not found"))
  }
  
  gene_expr <- expres_norm[gene_to_plot, ]
  
  # Recortar valores extremos para contraste visual
  max_cap <- quantile(gene_expr[gene_expr > 0], 0.99, na.rm = TRUE)
  gene_expr[gene_expr > max_cap] <- max_cap
  
  # Reescalar entre 0 y 1
  gene_expr_scaled <- scales::rescale(gene_expr, to = c(0, 1))
  
  # Asegurar que meta$X esté presente y sea índice correcto
  if (!"X" %in% colnames(meta)) {
    meta$X <- colnames(exprs_norm)
  }
  
  meta$gene_expr <- gene_expr_scaled[meta$X]
  
  # Ordenar metadata para que los puntos con baja expresión (gris) se dibujen primero
  meta <- meta[order(meta$gene_expr, decreasing = FALSE), ]
  
  # Plot
  p <- ggplot(meta, aes(x = UMAP_1, y = UMAP_2, color = gene_expr)) +
    geom_point(color = "black", size = 0.9, alpha = 0.6) +
    geom_point(size = 2) +
    #scale_color_gradientn(colors = c("lightgrey", "red")) +
    scale_color_viridis_c(direction = 1, na.value = "white") +
    labs(
      title = paste("Expression of", gene_to_plot),
      color = "Scaled expression"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank()
    ) +
    coord_fixed()
  
  # Si se seleccionó modo separado, usar facet_wrap
  if (split_mode == "Separated") {
    p <- p + facet_wrap(~condition)
  }
  
  return(p)
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
    "CTRL" = "#9932CC", 
    "SSc" = "#228B22",
    "Tn"= "#a9dfbb", 
    "Tcm"= "#6cd0c4", 
    "Treg"= "#fc6daf", 
    "Tem"= "#3faccd",
    "Temra"= "#ff0000", 
    "Th17"= '#d7ff55', 
    "Th2"= '#FF9F55', 
    "Thf"= '#559FFF',
    "Th1_GZMK"= '#55FF9F', 
    "Th1_GZMB"= '#FFC255',
    "Th22"='#55FF9F',
    "Th1"='#FFC255',
    "nTreg" = "#6A4C93",
    "iTreg" = "#008B8B",
    "eTreg" ="#B22222"
  )
  
  if (group_var == "cell_type") {
    df[[group_var]] <- factor(df[[group_var]], levels = c("Tn", "Tcm", "Tem", "Temra", "Treg"))
  } else if (group_var == "Th") {
      df[[group_var]] <- factor(df[[group_var]], levels = c("Tn", "Th1", "Th2", "Th17","Th22", "Thf","Tem", "Temra", "Treg"))
  } else if (group_var == "clusters") {
    df[[group_var]] <- factor(df[[group_var]], levels = as.character(0:11))
  } else if (group_var == "condition") {
    df[[group_var]] <- factor(df[[group_var]], levels = c("CTRL", "SSc"))
  } else if (group_var == "ct_Treg") {
    df[[group_var]] <- factor(df[[group_var]], levels = c("Tn", "Tcm","Tem", "Temra", "nTreg","iTreg","eTreg"))
  } else if (group_var == "ct_sub") {
    df[[group_var]] <- factor(df[[group_var]], levels = c("Tn", "Th1", "Th2", "Th17","Th22", "Thf","Tem", "Temra", "nTreg","iTreg","eTreg"))
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
    labs(title = paste("Mean expression of", gene),
         x = "",
         y = "",
         fill = "") +
    scale_fill_manual(values = custom_colors, na.translate = FALSE) +
    theme_bw() +
    theme(legend.position = 'none',
          axis.text.x = element_text(size = 12, face = "bold"),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_blank(),
          title = element_text(size = 20),
          plot.margin = unit(c(1, 1, 1, 1), "mm"))
}

# Function 2: Grouped by variable (e.g., R09 or cell_type)
plot_boxplot_expr_grouped <- function(gene, group_var, exprs_norm, meta) {
  meta$clusters <- as.character(meta$R09)
  gene_expr <- as.numeric(Matrix::as.matrix(exprs_norm[gene, , drop = FALSE]))
  names(gene_expr) <- colnames(exprs_norm)
  
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
    "CTRL" = "#9932CC", 
    "SSc" = "#228B22",
    "Tn"= "#a9dfbb", 
    "Tcm"= "#6cd0c4", 
    "Treg"= "#fc6daf", 
    "Tem"= "#3faccd",
    "Temra"= "#ff0000", 
    "Th17"= '#d7ff55', 
    "Th2"= '#FF9F55', 
    "Thf"= '#559FFF',
    "Th1_GZMK"= '#55FF9F', 
    "Th1_GZMB"= '#FFC255',
    "Th22"='#55FF9F',
    "Th1"='#FFC255',
    "nTreg" = "#6A4C93",
    "iTreg" = "#008B8B",
    "eTreg" ="#B22222"
  )
  
  df <- meta %>%
    mutate(expression = gene_expr[colnames(exprs_norm)]) %>%
    select(Sample, condition, clusters, cell_type, Th, ct_Treg, ct_sub, expression)
  
  group_vars <- c("Sample", group_var)
  if (group_var != "condition") group_vars <- c(group_vars, "condition")
  
  df_summary <- df %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(mean_expr = mean(expression, na.rm = TRUE), .groups = "drop")
  
  df_summary$facet_group <- factor(df_summary[[group_var]],
                                   levels = if (group_var == "cell_type") {
                                     c("Tn", "Tcm", "Tem", "Temra", "Treg")
                                     #annot <- "Cell type annotation"
                                   } else if (group_var == "Th") {
                                     c("Tn", "Th1", "Th2", "Th17","Th22", "Thf","Tem", "Temra", "Treg")
                                     #annot <- "Cell type T helper annotation"
                                   } else if (group_var == "ct_Treg") {
                                     c("Tn", "Tcm","Tem", "Temra", "nTreg","iTreg","eTreg")
                                     #annot <- "Cell type Treg annotation"
                                   } else if (group_var == "ct_sub") {
                                     c("Tn", "Th1", "Th2", "Th17","Th22", "Thf","Tem", "Temra", "nTreg","iTreg","eTreg")
                                     #annot <- "Cell type sublcustering annotation"
                                   } else if (group_var == "clusters"){
                                     as.character(0:11)
                                     #annot <- "Clusters"
                                   } else {
                                     c("CTRL","SSc")
                                     #annot <- "SSc vs CTRL"
                                   })

  
  p <- ggplot(df_summary, aes(x = condition, y = mean_expr, fill = .data[[group_var]])) +
    geom_boxplot(width = 0.5, alpha = 0.75, outlier.shape = NA) +
    geom_point(alpha = 0.9, shape = 21, colour = "black", stroke = 0.3,
               size = 2.5, position = position_jitter(width = 0.15)) +
    labs(title = paste0("Mean expression of ", gene),
         x = group_var, 
         y = "", 
         fill = group_var) +
    scale_fill_manual(values = custom_colors, na.translate = FALSE) +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 12, face = "bold", angle = 45, hjust = 1),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_blank(),
          title = element_text(size = 18),
          plot.margin = unit(c(1, 1, 1, 1), "mm"))
  
  if (group_var != "condition") {
    p <- p + facet_wrap(~facet_group)
  }
  
  p
}

# Heatmap ####
#' Infer metadata filtering rules from DEA dataset name
#' @param dea_name Character string with DEA dataset name
#' @return List with group_var and group_level
infer_heatmap_filter <- function(dea_name) {
  
  nm <- as.character(dea_name)
  nm_clean <- gsub("\\s+", "", nm)  # remove spaces
  
  # Global / no filtering
  if (tolower(nm_clean) %in% c("global", "all", "overall")) {
    return(list(group_var = NULL, group_level = NULL))
  }
  
  # Cluster pattern (e.g., Cluster0, clusters_3, cluster-5)
  if (grepl("cluster", nm_clean, ignore.case = TRUE)) {
    cl_num <- sub(".*?([0-9]+)$", "\\1", nm_clean)
    return(list(group_var = "clusters",
                group_level = as.character(cl_num)))
  }
  
  # Condition
  if (nm_clean %in% c("CTRL", "SSc")) {
    return(list(group_var = "condition",
                group_level = nm_clean))
  }
  
  # Cell type annotation
  if (nm_clean %in% c("Tn","Tcm","Tem","Temra","Treg")) {
    return(list(group_var = "cell_type",
                group_level = nm_clean))
  }
  
  # T helper annotation
  if (nm_clean %in% c("Th1","Th2","Th17","Th22","Thf")) {
    return(list(group_var = "Th",
                group_level = nm_clean))
  }
  
  # Treg subtype annotation
  if (nm_clean %in% c("nTreg","iTreg","eTreg")) {
    return(list(group_var = "ct_Treg",
                group_level = nm_clean))
  }
  
  # Fallback: no filtering
  return(list(group_var = NULL, group_level = NULL))
}

#' Scale rows using z-score normalization and clip extreme values
#' @param x A numeric matrix (genes x samples)
#' @return Scaled matrix with values capped between -2 and 2
scale_rows <- function(x) {
  x_scaled <- t(scale(t(x)))
  x_scaled[x_scaled > 2] <- 2
  x_scaled[x_scaled < -2] <- -2
  return(x_scaled)
}

#' Compute average expression per sample
#' @param expr_mat Gene expression matrix (genes x cells)
#' @param meta_samples Vector of sample IDs per cell (must align with columns of expr_mat)
#' @param top_genes Vector of gene names to retain
#' @return Matrix (genes x samples) with averaged expression
average_by_sample <- function(expr_mat, meta_samples, top_genes) {
  samples_unique <- unique(meta_samples)
  exp_ave_mat <- matrix(NA, nrow = length(top_genes), ncol = length(samples_unique))
  rownames(exp_ave_mat) <- top_genes
  colnames(exp_ave_mat) <- samples_unique
  
  for (i in seq_along(samples_unique)) {
    sample_id <- samples_unique[i]
    cells <- which(meta_samples == sample_id)
    if (length(cells) > 1) {
      exp_ave_mat[, i] <- rowMeans(expr_mat[, cells, drop = FALSE])
    } else {
      exp_ave_mat[, i] <- expr_mat[, cells]
    }
  }
  return(exp_ave_mat)
}

#' Rename sample columns using a named vector
#' @param matrix A matrix with sample column names
#' @param name_map Named vector: names are original sample names, values are new names
#' @return Matrix with updated column names
rename_samples <- function(matrix, name_map) {
  colnames(matrix) <- ifelse(
    colnames(matrix) %in% names(name_map),
    name_map[colnames(matrix)],
    colnames(matrix)
  )
  return(matrix)
}

#' Generate a heatmap-ready matrix from expression and metadata
#' @param exprs_norm Normalized gene expression matrix (genes x cells)
#' @param meta Metadata containing Sample column
#' @param top_genes Vector of genes to plot
#' @param sample_names Dataframe with columns og and new (original and new sample names)
#' @return Scaled and renamed heatmap matrix (genes x samples)
#' Generate heatmap-ready matrix with metadata filtering
#' @param exprs_norm Normalized gene expression matrix (genes x cells)
#' @param meta Filtered metadata
#' @param top_genes Vector of genes to plot
#' @param sample_names Dataframe with original and new sample names
#' @return Scaled and renamed matrix (genes x samples)
generate_heatmap_matrix <- function(exprs_norm,
                                    meta,
                                    top_genes,
                                    sample_names) {
  
  # Identify cell IDs to retain
  if ("X" %in% colnames(meta)) {
    cells_to_keep <- meta$X
  } else {
    cells_to_keep <- rownames(meta)
  }
  
  # Ensure intersection with expression matrix
  cells_to_keep <- intersect(cells_to_keep,
                             colnames(exprs_norm))
  
  # Subset expression matrix
  expr_sub <- Matrix::as.matrix(
    exprs_norm[top_genes, cells_to_keep, drop = FALSE]
  )
  
  # Extract sample labels
  meta_samples <- as.character(meta$Sample)
  
  # Compute average per sample
  exp_ave_mat <- average_by_sample(expr_sub,
                                   meta_samples,
                                   top_genes)
  
  # Scale rows (z-score clipping)
  exp_ave_scale <- scale_rows(exp_ave_mat)
  
  # Rename samples
  name_map <- setNames(sample_names$new,
                       sample_names$og)
  
  exp_ave_scale <- rename_samples(exp_ave_scale,
                                  name_map)
  
  return(exp_ave_scale)
}
#
# UMAP ####
# Function to generate UMAP plots
plot_umap <- function(meta, color_by = "clusters", split_mode = "Combined") {

  color_pal <- c(
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
    "CTRL" = "#9932CC", 
    "SSc" = "#228B22",
    "Tn"= "#a9dfbb", 
    "Tcm"= "#6cd0c4", 
    "Treg"= "#fc6daf", 
    "Tem"= "#3faccd",
    "Temra"= "#ff0000", 
    "Th17"= '#d7ff55', 
    "Th2"= '#FF9F55', 
    "Thf"= '#559FFF',
    "Th1_GZMK"= '#55FF9F', 
    "Th1_GZMB"= '#FFC255',
    "Th22"='#55FF9F',
    "Th1"='#FFC255',
    "nTreg" = "#6A4C93",
    "iTreg" = "#008B8B",
    "eTreg" ="#B22222"
  )
  
  meta$clusters <- as.character(meta$R09)
  meta$umap_density_condition <- as.numeric(meta$umap_density_condition)
  meta$cc_aa_size <- as.numeric(meta$cc_aa_size)
  umap_data <- meta %>% mutate(var_to_plot = .[[color_by]])
  is_cont <- is.numeric(umap_data$var_to_plot)
  
  if(color_by == "clusters"){
    annot <- "Clusters"
  } else if (color_by == "cell_type") {
    annot <- "Cell type annotation"
  } else if (color_by == "Th") {
    annot <- "T helper annotation"
  } else if (color_by == "ct_Treg") {
    annot <- "Treg annotation"
  } else if (color_by == "ct_sub" ) {
    annot <- "Subclustering annotatition"
  } else if (color_by == "condition") {
    annot <- "SSc vs CTRL"
  } else if (color_by == "Sample") {
    annot <- "samples"
  } else if (color_by == "percent_mito") {
    annot <- "Percentage of Mitochondrial content"
  } else if (color_by == "umap_density_condition") {
    annot <- "Cell density"
  } else if (color_by == "cc_aa_size") {
    annot <- "Clonotype expansion"
  } 
  
  if (is_cont) {
    umap_data <- umap_data %>%
      dplyr::arrange(is.na(var_to_plot), var_to_plot)
  }
  
  base_plot <- ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, color = var_to_plot)) +
    geom_point(color = "black", size = 0.9, alpha = 0.8) +
    geom_point(aes(color = var_to_plot), size = 0.7) +
    labs(title = paste("UMAP colored by", annot), color = "") +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.text = element_blank(), 
      axis.ticks = element_blank(),
      axis.title = element_text(size = 10),
      panel.grid = element_blank(), 
      panel.border = element_rect(color = "black", size = 1),
      legend.position = "right", 
      legend.text = element_text(size = 18)
    )+
    guides(color = guide_legend(override.aes = list(size = 4)))
  
  vars_with_custom_palette <- c("clusters", "condition", "cell_type","Th","ct_Treg","ct_sub")
  
  if (is_cont) {
    base_plot <- base_plot + scale_color_viridis_b(option = "plasma", direction = 1, na.value = "grey80")
  } else if (color_by %in% vars_with_custom_palette) {
    base_plot <- base_plot + scale_color_manual(values = color_pal, na.value = "grey80")
  } else {
    levs <- sort(unique(umap_data$var_to_plot))
    n_colors <- length(levs)
    palette <- suppressWarnings(RColorBrewer::brewer.pal(min(max(3, n_colors), 8), "Set2"))
    if (n_colors > 8) {
      palette <- colorRampPalette(palette)(n_colors)
    }
    base_plot <- base_plot + scale_color_manual(values = setNames(palette, levs), na.value = "grey80")
  }
  
  base_plot <- base_plot + coord_fixed()
  
  if (split_mode == "Separated") {
    base_plot <- base_plot + facet_wrap(~condition)
  }
  
  return(base_plot)
}
#
# Dotplot ####
plot_cluster_dotplot <- function(expr_mat,
                                 metadata,
                                 marker_genes,
                                 cluster_col,
                                 split_mode = "Combined") {
  
  stopifnot(inherits(expr_mat, "dgCMatrix"))
  stopifnot(cluster_col %in% colnames(metadata))
  
  marker_genes <- intersect(marker_genes, rownames(expr_mat))
  if (length(marker_genes) == 0) {
    stop("No marker genes found in expression matrix.")
  }
  
  expr_sub <- expr_mat[marker_genes, , drop = FALSE]
  
  clusters  <- as.character(metadata[[cluster_col]])
  
  # ------------------------------------------------------------
  # CASE 1 — Combined mode
  # ------------------------------------------------------------
  if (split_mode == "Combined") {
    
    group_ids <- clusters
    group_df  <- data.frame(group = group_ids)
    
  } else {
    
    # ------------------------------------------------------------
    # CASE 2 — Separated mode (cluster + condition)
    # ------------------------------------------------------------
    
    stopifnot("condition" %in% colnames(metadata))
    
    conditions <- as.character(metadata$condition)
    group_ids  <- paste(clusters, conditions, sep = "__")
    
    group_df <- data.frame(
      group = group_ids,
      cluster = clusters,
      condition = conditions
    )
  }
  
  group_levels <- unique(group_ids)
  
  result_list <- lapply(group_levels, function(gr) {
    
    cells_idx <- which(group_ids == gr)
    if (length(cells_idx) == 0) return(NULL)
    
    mat_cl <- expr_sub[, cells_idx, drop = FALSE]
    
    avg_expr <- Matrix::rowMeans(mat_cl)
    pct_expr <- Matrix::rowSums(mat_cl > 0) / length(cells_idx) * 100
    
    data.frame(
      group = gr,
      gene = names(avg_expr),
      avg_expr = as.numeric(avg_expr),
      pct_expr = as.numeric(pct_expr),
      stringsAsFactors = FALSE
    )
  })
  
  dot_data <- do.call(rbind, result_list)
  
  # ------------------------------------------------------------
  # Extract cluster + condition back if separated
  # ------------------------------------------------------------
  
  if (split_mode == "Separated") {
    
    split_info <- do.call(rbind, strsplit(dot_data$group, "__"))
    dot_data[[cluster_col]] <- split_info[,1]
    dot_data$condition <- split_info[,2]
    
  } else {
    
    dot_data[[cluster_col]] <- dot_data$group
  }
  
  # ------------------------------------------------------------
  # Z-score per gene
  # ------------------------------------------------------------
  
  dot_data <- dot_data %>%
    dplyr::group_by(gene) %>%
    dplyr::mutate(avg_expr_scaled = as.numeric(scale(avg_expr))) %>%
    dplyr::ungroup()
  
  dot_data$avg_expr_scaled[is.na(dot_data$avg_expr_scaled)] <- 0
  
  # ------------------------------------------------------------
  # Factor ordering
  # ------------------------------------------------------------
  
  if (cluster_col == "cell_type") {
    dot_data[[cluster_col]] <- factor(dot_data[[cluster_col]],
                                      levels = c("Tn","Tcm","Tem","Temra","Treg"))
    annot <- "cell type annotation"
  } else if (cluster_col == "Th") {
    dot_data[[cluster_col]] <- factor(dot_data[[cluster_col]],
                                      levels = c("Tn","Th1","Th2","Th17","Th22","Thf","Tem","Temra","Treg"))
    annot <- "T helper annotation"
  } else if (cluster_col == "clusters") {
    dot_data[[cluster_col]] <- factor(dot_data[[cluster_col]],
                                      levels = as.character(0:11))
    annot <- "Clusters"
  } else if (cluster_col == "condition") {
    dot_data[[cluster_col]] <- factor(dot_data[[cluster_col]],
                                      levels = c("CTRL","SSc"))
    annot <- "SSc vs CTRL"
  } else if (cluster_col == "ct_Treg") {
    dot_data[[cluster_col]] <- factor(dot_data[[cluster_col]], 
                                      levels = c("Tn", "Tcm","Tem", "Temra", "nTreg","iTreg","eTreg"))
    annot <- "Treg annotation"
  } else if (cluster_col == "ct_sub") {
    dot_data[[cluster_col]] <- factor(dot_data[[cluster_col]], 
                                      levels = c("Tn", "Th1", "Th2", "Th17","Th22", "Thf","Tem", "Temra", "nTreg","iTreg","eTreg"))
    annot <- "Subclustering annotation"
  } 
  
  # ------------------------------------------------------------
  # Plot
  # ------------------------------------------------------------
  
  p <- ggplot(dot_data,
              aes(x = .data[[cluster_col]], y = gene)) +
    geom_point(aes(size = pct_expr,
                   fill = avg_expr_scaled),
               shape = 21,
               colour = "black",
               stroke = 0.3) +
    scale_fill_viridis_c(direction = 1, na.value = "white") +
    scale_size(range = c(1, 7)) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12,face = "bold"),
      axis.text.y = element_text(size = 12,face = "bold"),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    ) +
    labs(
      title = paste0("Dotplot grouped by ", annot),
      fill  = "Scaled avg expr",
      size  = "% expressing"
    ) +
    coord_flip()
  
  if (split_mode == "Separated") {
    p <- p + facet_wrap(~condition)
  }
  
  return(p)
}
#
# Demographic ----
generate_demographic_summary <- function(demo_df) {
  
  n_total <- nrow(demo_df)
  n_patients <- n_total/2
  
  # Sex
  n_female <- sum(demo_df$Sex == "Female", na.rm = TRUE)
  pct_female <- round(n_female / n_total * 100, 1)
  
  # Age
  age_mean_ssc <- round(mean((demo_df %>% filter(condition == "SSc"))$Age, na.rm = TRUE), 1)
  age_sd_ssc   <- round(sd((demo_df %>% filter(condition == "SSc"))$Age, na.rm = TRUE), 1)
  
  age_mean_ctrl <- round(mean((demo_df %>% filter(condition == "CTRL"))$Age, na.rm = TRUE), 1)
  age_sd_ctrl   <- round(sd((demo_df %>% filter(condition == "CTRL"))$Age, na.rm = TRUE), 1)
  
  summary_list <- list(
    data.frame(
      Variable = "N patients",
      Value = as.character(n_patients)
    ),
    data.frame(
      Variable = "N controls",
      Value = as.character(n_patients)
    ),
    data.frame(
      Variable = "Female (%)",
      Value = paste0(n_female, " (", pct_female, "%)")
    ),
    data.frame(
      Variable = "SSc Age (mean ± SD)",
      Value = paste0(age_mean_ssc, " ± ", age_sd_ssc)
    ),
    data.frame(
      Variable = "CTRL Age (mean ± SD)",
      Value = paste0(age_mean_ctrl, " ± ", age_sd_ctrl)
    )
  )
  
  clinical_cols <- setdiff(colnames(demo_df),
                           c("Sample", "Sex", "Age"))
  
  clinical_summary <- lapply(clinical_cols, function(col) {
    
    demo_df <- demo_df %>% filter(condition == "SSc")
    
    values <- demo_df[[col]]
    
    n_pos <- sum(toupper(values) == "YES", na.rm = TRUE)
    n_valid <- sum(toupper(values) %in% c("YES", "NO"))
    
    pct <- ifelse(n_valid > 0,
                  round(n_pos / n_valid * 100, 1),
                  NA)
    
    data.frame(
      Variable = col,
      Value = ifelse(is.na(pct),
                     paste0(n_pos),
                     paste0(n_pos, " (", pct, "%)"))
    )
  })
  
  summary_df <- dplyr::bind_rows(summary_list,
                                 dplyr::bind_rows(clinical_summary))
  
  return(summary_df)
}

plot_demographic_barplot <- function(demo_df) {
  
  demo_df <- demo_df %>% filter(condition == "SSc")
  
  clinical_cols <- setdiff(colnames(demo_df),
                           c("Sample", "Sex", "Age"))
  
  df_plot <- lapply(clinical_cols, function(col) {
    
    values <- toupper(demo_df[[col]])
    
    n_valid <- sum(values %in% c("YES", "NO"))
    n_pos   <- sum(values == "YES", na.rm = TRUE)
    
    pct <- ifelse(n_valid > 0,
                  n_pos / n_valid * 100,
                  NA)
    
    data.frame(
      Variable = col,
      Percentage = pct
    )
    
  }) |> dplyr::bind_rows()
  
  df_plot <- df_plot |>
    dplyr::filter(!is.na(Percentage)) |>
    dplyr::arrange(Percentage)
  
  ggplot(df_plot,
         aes(x = Percentage,
             y = reorder(Variable, Percentage))) +
    geom_col(fill = "#2c7fb8") +
    labs(x = "Percentage of patients (%)",
         y = "",
         title = "Clinical characteristics of the SSc cohort") +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid = element_blank()
    )
}

#
<<<<<<< HEAD
# Libraries, functions and data ####

library("shiny")
library("devtools")
library("tidyverse")
library("ggpubr")
library("ggrepel")
library("matrixStats")
library("RColorBrewer")
library("reshape2")
library("DT")
library("ggpubr")
library("DelayedArray")
library("shinyjs")
library("cowplot")
library("pheatmap")
library("ggdendro")
library("ggsignif")
library("bslib")
library("rstatix")
library("enrichR")
library("Seurat")
set.seed(42)

source("functions/functionsShiny.R")
load("data/SScCD4_data_v2.RData")

=======
# Functions and data
data_lib <- c("devtools", "tidyverse", "ggpubr", "ggrepel", "enrichR",
              "RColorBrewer", "reshape2","shiny","DT","ggpubr",
              "bslib","shinyjs","ggsignif","ggdendro","cowplot","pheatmap","rstatix", "Seurat")
suppressPackageStartupMessages({
  invisible(lapply(data_lib, require, character.only = TRUE))
})
set.seed(42)

source("functions/functionsShiny.R")
source("functions/heatmapFunction.R")

load("data/SScCD4_data.RData")
load("data/dea_list.RData")
load("data/gseaTF_list.RData")
>>>>>>> f40c6eb162999aea48b4d3f321b5e2bce1953d15

# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # theme
  theme = bslib::bs_theme(version = 4, bootswatch = "litera"),
  
  # Logos
  column(
    12,
    tags$div(
      style = "display: flex; justify-content: center; align-items: center; gap: 40px; padding: 10px;",
<<<<<<< HEAD
      tags$img(src = "MainPage/project_logo.png", style = "height:170px;"),
=======
      tags$img(src = "MainPage/project_logo.png", style = "height:150px;"),
>>>>>>> f40c6eb162999aea48b4d3f321b5e2bce1953d15
      tags$img(src = "MainPage/ipblnlogo.jpg", style = "height:60px;"),
      tags$img(src = "MainPage/BCH-logo.png", style = "height:70px;"),
      tags$img(src = "MainPage/Broad-logo.jpeg", style = "height:40px;"),
      tags$img(src = "MainPage/MGA_logo_v6.png", style = "height:100px;")
    )
  ),

    # Application title
<<<<<<< HEAD
    #titlePanel("SSc CD4 project"),
=======
    titlePanel("SSc CD4 project"),
>>>>>>> f40c6eb162999aea48b4d3f321b5e2bce1953d15
  
  #ui Table####

  
  tabsetPanel(
    tabPanel("Differential Expression analysis",

             sidebarLayout(
               sidebarPanel(
                 fluidRow(
                   column(12,
<<<<<<< HEAD
                          h4("Select Differential expression analysis SSc vs CTRL"),
                          selectInput("dea_source", "",
=======
                          h4("Select Differential expression analysis"),
                          selectInput("dea_source", "Select DEA dataset:",
>>>>>>> f40c6eb162999aea48b4d3f321b5e2bce1953d15
                                      choices = names(dea_list),
                                      selected = "pseudobulk")
                   )
                 ),
                 width = 12,
                 actionButton(
                   inputId = "runAnalysis",
                   label = "Show table",
                   class = "btn-block small-btn",
                   style = "width: 10ß0%; padding: 4px 4px; font-size: 16px;"
                 ),
               ),
               mainPanel(width = 12,
                 tabPanel("DEA table", DTOutput("deaTable")),
                 )
               )
             ),

    
    # ui Volcano ####
    
    tabPanel("Volcano Plot",
             fluidRow(
               column(12,
                      tags$div(
                        style = "padding: 0 15px;",
<<<<<<< HEAD
                        h4("Volcano plot (Global DEA only)")
=======
                        selectInput("dea_source_volcano", "Select DEA dataset:",
                                    choices = names(dea_list),
                                    selected = "pseudobulk",
                                    width = "100%")
>>>>>>> f40c6eb162999aea48b4d3f321b5e2bce1953d15
                      )
               )
             ),
             width = 12,
             sidebarLayout(
               sidebarPanel(
                 selectizeInput("genes_heatmap_volcano", 
                                "Genes to label:",
                                choices = NULL,
                                multiple = TRUE,
                                selected = c("STAT1", "STAT3"),
                                options = list(maxItems = 50)),
<<<<<<< HEAD
                 sliderInput("pval_cutoff_volcano", "FDR threshold:", min = 0, max = 1, value = 0.05, step = 0.01),
=======
                 sliderInput("pval_cutoff_volcano", "FDR threshold:", min = 0, max = 0.2, value = 0.05, step = 0.01),
>>>>>>> f40c6eb162999aea48b4d3f321b5e2bce1953d15
                 actionButton("runVolcanoOnly", "Run Volcano Plot")
               ),
               mainPanel(
                 plotOutput("volcanoPlot", height = "500px")
               )
             )
    ),
  
    # ui Heatmap ####  
    tabPanel("Heatmap",
             fluidPage(
               fluidRow(
                 column(6,
                        h4("Select Differential expression analysis"),
                        selectInput("dea_source_heatmap", "Select DEA dataset:",
                                    choices = names(dea_list),
<<<<<<< HEAD
                                    selected = "Global",
                                    width = "auto"),
=======
                                    selected = "pseudobulk",
                                    width = "100%"),
>>>>>>> f40c6eb162999aea48b4d3f321b5e2bce1953d15
                        selectizeInput("genes_heatmap_only", 
                                       "Genes to plot:",
                                       choices = NULL,
                                       multiple = TRUE,
                                       selected = NULL,
                                       options = list(maxItems = 50),
<<<<<<< HEAD
                                       width = "auto"),
                        actionButton("runHeatmap", "Run Heatmap",
                                     style = "width: 10ß0%; padding: 6px; font-size: 16px;")
=======
                                       width = "100%"),
                        actionButton("runHeatmap", "Run Heatmap",
                                     style = "width: 100%; padding: 6px; font-size: 16px;")
>>>>>>> f40c6eb162999aea48b4d3f321b5e2bce1953d15
                 ),
                 column(8,
                        div(
                          style = "overflow-x: auto; overflow-y: auto; max-height: 1000px;",
                          imageOutput("heatmapPlot", width = "100%", height = "auto")
                        )
                 )
               )
             )
    ),
    
    #ui UMAP ####
    
    
    tabPanel("UMAP",
             sidebarLayout(
<<<<<<< HEAD
               sidebarPanel(width = 12,
                 selectInput("umap_color_by", "Color by:",
                             choices = c("clusters", "cell_type", "condition", "Sample",
                                         "percent_mito"), #"cellType_Th",
                             selected = "clusters")
               ),
               mainPanel(width = 12,
=======
               sidebarPanel(
                 selectInput("umap_color_by", "Color by:",
                             choices = c("clusters", "cell_type", "cellType_Th","condition", "Sample",
                                         "percent_mito"),
                             selected = "clusters")
               ),
               mainPanel(
>>>>>>> f40c6eb162999aea48b4d3f321b5e2bce1953d15
                 plotOutput("umapPlot", height = "600px")
               )
             )
    ),
    
    
<<<<<<< HEAD
    #ui UMAP expression#####
    
    tabPanel("Gene expression in UMAP",
             sidebarLayout(
               sidebarPanel(width = 12,
=======
    #ui MAP expression#####
    
    tabPanel("Gene expression in UMAP",
             sidebarLayout(
               sidebarPanel(
>>>>>>> f40c6eb162999aea48b4d3f321b5e2bce1953d15
                 selectizeInput(
                   options = list(placeholder = 'Type a gene name...'),
                   inputId = "gene_umap",
                   label = "Select gene:",
<<<<<<< HEAD
                   selected = "CD4",
                   choices = NULL,
                   multiple = FALSE
                 )
               ),
               mainPanel(width = 12,
=======
                   choices = NULL,
                   multiple = FALSE,
                   selected = "CD4"
                 )
               ),
               mainPanel(
>>>>>>> f40c6eb162999aea48b4d3f321b5e2bce1953d15
                 plotOutput("umapPlot_expression", height = "600px")
               )
             )
    ),
    

    #ui Box plot#####
    
    tabPanel("Gene Expression Boxplot",
             fluidRow(
               column(6,
<<<<<<< HEAD
                      selectizeInput("gene_boxplot", "Select gene:", 
                                     choices = NULL,  
                                     selected = "CD4",
                                     options = list(maxOptions = 5000))
=======
                      selectInput("gene_boxplot", "Select gene:", 
                                  choices = rownames(exprs_norm), 
                                  selected = "CD4")
>>>>>>> f40c6eb162999aea48b4d3f321b5e2bce1953d15
               ),
               column(6,
                      selectInput("group_by", "Group by:",
                                  choices = c("condition", "clusters", "cell_type"),
                                  selected = "condition")
               )
             ),
             fluidRow(
               column(12, h4("Per-sample expression boxplot")),
               column(12, plotOutput("boxplot_expr_individual",height = "200px"))
             ),
             fluidRow(
               column(12, h4("Grouped boxplot (per group + condition)")),
               column(12, plotOutput("boxplot_expr_grouped",height = "600px"))
             )
    ),
    
    # ui GSEA ####
    tabPanel("GSEA",
             fluidRow(
               column(4,
                      radioButtons("gsea_category", "Gene set category:",
<<<<<<< HEAD
                                   choices = c("Hallmark" = "H", "Reactome" = "Reactome"),
=======
                                   choices = c("Hallmark" = "H", "C2" = "C2"),
>>>>>>> f40c6eb162999aea48b4d3f321b5e2bce1953d15
                                   selected = "H", inline = TRUE)
               ),
               column(4,
                      selectInput("gsea_dataset", "Select GSEA dataset:",
<<<<<<< HEAD
                                  choices = names(gseaTF_list), 
                                  selected = if ("Global_Hallmark" %in% names(gseaTF_list)) "Global_Hallmark" else names(gseaTF_list)[1])
               ),
               column(4,
                      numericInput("fdr_gsea", "FDR cutoff:", value = 0.1, step = 0.01),
=======
                                  choices = names(gseaTF_list), selected = NULL)
               ),
               column(4,
                      numericInput("fdr_gsea", "FDR cutoff:", value = 0.05, step = 0.01),
>>>>>>> f40c6eb162999aea48b4d3f321b5e2bce1953d15
                      numericInput("top_n_gsea", "Top pathways per direction:", value = 10, step = 1)
               )
             ),
             
             tabsetPanel(
               tabPanel("Table",
                        DTOutput("gsea_table")
               ),
               tabPanel("Plot",
                        plotOutput("gsea_plot", height = "600px")
               )
             )
    )
    
    
    ),
  

  
  
)


# server ####
# Define server logic required to draw a histogram
server <- function(input, output,session) {
  
  output$About<-renderText({
    HTML('<style> .div-1 {background-color: #EBEBEB;} </style>
    <div-1 style="text-align: center;" class="container">
    <h4>About this website</h4>
    <p> </p>
    <p>Analysis developed by Gonzalo Villanueva-Martin.</p>
    <p>This website was developed by Gonzalo Villanueva-Martin.</p>
    <p>For any questions or comments please contact Gonzalo Villanueva-Martin , 
    </p>
    </div>
    ')
  })
  
<<<<<<< HEAD


  
  #server Table#####
  
  dea_data <- reactive({
    dea <- dea_list[[input$dea_source]]
    req(dea)
    
    # ---- Harmonize column names depending on format ---- #
    
    if (all(c("Estimate", "adj.P.Val") %in% colnames(dea))) {
      # Global-style (already OK)
      
      dea_std <- dea %>%
        mutate(
          logFC = Estimate,
          padj  = adj.P.Val
        )
      
    } else if (all(c("logFC", "padj") %in% colnames(dea))) {
      # Cluster-style (new format)
      
      dea_std <- dea
      
    } else {
      validate(
        need(FALSE, "DEA table has unknown column format")
      )
    }
    
    # ---- Add significance consistently ---- #
    
    dea_std <- dea_std %>%
      mutate(
        significance = case_when(
          padj < 0.05 & logFC > 0 ~ "Upregulated",
          padj < 0.05 & logFC < 0 ~ "Downregulated",
          TRUE ~ "Not_significant"
        )
      )
    
    return(dea_std)
  })

  
  # Observer for table
  observeEvent(input$runAnalysis, {
    output$deaTable <- renderDT({
      datatable(dea_data(), 
                options = list(pageLength = 20, scrollX = TRUE, autoWidth = TRUE), 
=======
  dea_data <- reactive({
    dea <- dea_list[[input$dea_source]]
    
    # Asegurar columnas mínimas
    req(all(c("gene", "Estimate", "adj.P.Val") %in% colnames(dea)))
    
    # Asegurar columna de "significance" si no existe
    if (!"significance" %in% colnames(dea)) {
      dea <- dea %>%
        mutate(
          significance = case_when(
            adj.P.Val < 0.05 & Estimate > 0 ~ "Upregulated",
            adj.P.Val < 0.05 & Estimate < 0 ~ "Downregulated",
            TRUE ~ "Not_significant"
          )
        )
    }
    
    return(dea)
  })

  
  #server Table#####

  
  # Observador for table
  observeEvent(input$runAnalysis, {
    output$deaTable <- renderDT({
      datatable(dea_data(), 
                options = list(pageLength = 10, scrollX = TRUE), 
>>>>>>> f40c6eb162999aea48b4d3f321b5e2bce1953d15
                filter = 'top')
    })
  })
  
  
<<<<<<< HEAD
  #server Volcano #####
  
  observe({
    updateSelectizeInput(
      session,
      inputId = "genes_heatmap_volcano",
      choices = unique(dea_list[["Global"]]$gene),
=======
  #server Volcano and heatmap #####
  
  observe({
    req(dea_data())
    updateSelectizeInput(
      session,
      inputId = "genes_heatmap",
      choices = unique(dea_data()$gene),
      server = TRUE
    )
  })
  
  
  # Volcano plot only
  observe({
    req(input$dea_source_volcano)
    updateSelectizeInput(
      session,
      inputId = "genes_heatmap_volcano",
      choices = unique(dea_list[[input$dea_source_volcano]]$gene),
>>>>>>> f40c6eb162999aea48b4d3f321b5e2bce1953d15
      server = TRUE
    )
  })
  
  observeEvent(input$runVolcanoOnly, {
<<<<<<< HEAD
    
    volcano_data <- dea_list[["Global"]]
    
    volcano_data$significance <- case_when(
      volcano_data$adj.P.Val < input$pval_cutoff_volcano & volcano_data$Estimate > 0 ~ "Upregulated",
      volcano_data$adj.P.Val < input$pval_cutoff_volcano & volcano_data$Estimate < 0 ~ "Downregulated",
=======
    volcano_data <- dea_list[[input$dea_source_volcano]]
    
    volcano_data$significance <- case_when(
      volcano_data$adj.P.Val < 0.05 & volcano_data$Estimate > 0 ~ "Upregulated",
      volcano_data$adj.P.Val < 0.05 & volcano_data$Estimate < 0 ~ "Downregulated",
>>>>>>> f40c6eb162999aea48b4d3f321b5e2bce1953d15
      TRUE ~ "Not_significant"
    )
    
    volcano_data$delabel <- ifelse(
      volcano_data$gene %in% input$genes_heatmap_volcano,
      volcano_data$gene,
      NA
    )
    
<<<<<<< HEAD
    fdr_threshold <- -log10(max(volcano_data$p_value[volcano_data$adj.P.Val <= input$pval_cutoff_volcano], na.rm = TRUE))
    
    output$volcanoPlot <- renderPlot({
      
      volcano_data$highlight <- volcano_data$gene %in% input$genes_heatmap_volcano
      
      ggplot(volcano_data, aes(x = Estimate, y = -log10(p_value), label = delabel)) +
        
        # Normal points
        geom_point(alpha = 0.9, shape = 21, colour = "black", stroke = 0.3, size = 2.5,
                   aes(fill = significance)) +
        
        # Highlight points
        geom_point(
          data = filter(volcano_data, highlight),
          aes(text = delabel),
          color = "black", fill = "red",
          alpha = 0.9, shape = 21, stroke = 0.3, size = 6
        ) +
        
        # Theme
        theme_minimal() +
        theme(axis.title = element_text(size=16),
              legend.text = element_text(size=16),
              legend.title = element_blank())+
        labs(x = "Estimate", y = "-log10(p-value)") +
        
        # Labels
        geom_label_repel(data = filter(volcano_data, !is.na(delabel)),
                         color = "black", 
                         max.overlaps = 10, 
                         force = 5, 
                         size = 6,
                         show.legend = FALSE) +
        
        # Colors
        scale_fill_manual(values = c("Upregulated" = "#F77663",
                                     "Downregulated" = "#3D799E",
                                     "Not_significant" = "#f0f0f0")) +
        
        # Lines
=======
    fdr_threshold <- -log10(input$pval_cutoff_volcano)
    
    output$volcanoPlot <- renderPlot({
      ggplot(volcano_data, aes(x = Estimate, y = -log10(adj.P.Val), label = delabel)) +
        geom_point(alpha = 0.9, shape = 21, colour = "black", stroke = 0.3, size = 2.5,
                   aes(fill = significance)) +
        theme_minimal() +
        labs(x = "Estimate", y = "-log10(Adj p-value)") +
        geom_label_repel(data = filter(volcano_data, !is.na(delabel)),
                         color = "black", max.overlaps = 10, force = 5, size = 3,
                         show.legend = FALSE) +
        scale_fill_manual(values = c("Upregulated" = "#F77663",
                                     "Downregulated" = "#3D799E",
                                     "Not_significant" = "#f0f0f0")) +
>>>>>>> f40c6eb162999aea48b4d3f321b5e2bce1953d15
        geom_hline(yintercept = fdr_threshold, linetype = "dashed", color = "blue") +
        geom_vline(xintercept = 0, linetype = "dashed", color = "black")
    })
  })
  
  #server heatmap #####
<<<<<<< HEAD
  
  #exprs_matrix <- reactive({
    #readRDS("data/exprs_norm.rds")
  #})
  
=======
>>>>>>> f40c6eb162999aea48b4d3f321b5e2bce1953d15
  observe({
    req(input$dea_source_heatmap)
    updateSelectizeInput(
      session,
      inputId = "genes_heatmap_only",
      choices = unique(dea_list[[input$dea_source_heatmap]]$gene),
      server = TRUE
    )
  })
  
  observeEvent(input$runHeatmap, {
<<<<<<< HEAD
    req(input$dea_source_heatmap)
    
=======
>>>>>>> f40c6eb162999aea48b4d3f321b5e2bce1953d15
    dea_data <- dea_list[[input$dea_source_heatmap]]
    
    if (length(input$genes_heatmap_only) > 0) {
      top_genes <- input$genes_heatmap_only
    } else {
<<<<<<< HEAD
      DE_pos <- dea_data %>% filter(adj.P.Val < 0.05, Estimate > 0) %>%
        arrange(desc(Estimate)) %>% slice_head(n = 20)
      DE_neg <- dea_data %>% filter(adj.P.Val < 0.05, Estimate < 0) %>%
        arrange(Estimate) %>% slice_head(n = 20)
      top_genes <- unique(bind_rows(DE_pos, DE_neg)$gene)
    }
    
    #expr_mat <- exprs_matrix()
    exp_heat <- Matrix::as.matrix(exprs_norm[top_genes, , drop = FALSE])
    
    # Asegurar que las muestras estén alineadas con la metadata
    meta_samples <- as.character(meta$Sample)
    samples_unique <- unique(meta_samples)
    
    # Inicializar matriz de promedio por muestra
    exp_ave_mat <- matrix(NA, nrow = length(top_genes), ncol = length(samples_unique))
    rownames(exp_ave_mat) <- top_genes
    colnames(exp_ave_mat) <- samples_unique
    
    # Calcular la media por muestra
    for (i in seq_along(samples_unique)) {
      sample_id <- samples_unique[i]
      cells <- which(meta_samples == sample_id)
      if (length(cells) > 1) {
        exp_ave_mat[, i] <- rowMeans(exp_heat[, cells, drop = FALSE])
      } else {
        exp_ave_mat[, i] <- exp_heat[, cells]
      }
    }
    
    # Z-score row scaling
    scale_rows <- function(x) t(scale(t(x)))
    exp_ave_scale <- scale_rows(exp_ave_mat)
=======
      DE_pos <- dea_data %>%
        filter(adj.P.Val < 0.05, Estimate > 0) %>%
        arrange(desc(Estimate)) %>%
        slice_head(n = 20)
      
      DE_neg <- dea_data %>%
        filter(adj.P.Val < 0.05, Estimate < 0) %>%
        arrange(Estimate) %>%
        slice_head(n = 20)
      
      top_genes <- unique(bind_rows(DE_pos, DE_neg)$gene)
    }
    
    exp_heat <- exprs_norm[top_genes, ]
    exp_heat <- as.data.frame(t(exp_heat))
    exp_heat$sample <- as.character(meta$Sample)
    
    exp_ave <- aggregate(exp_heat[, 1:(ncol(exp_heat)-1)], list(exp_heat$sample), mean)
    colnames(exp_ave)[1] <- "sample"
    exp_ave <- as.data.frame(t(exp_ave))
    colnames(exp_ave) <- as.character(t(exp_ave[1, ]))
    exp_ave <- exp_ave[-1, ]
    row_names <- rownames(exp_ave)
    exp_ave <- mutate_all(exp_ave, function(x) as.numeric(as.character(x)))
    rownames(exp_ave) <- row_names
    exp_ave <- as.matrix(exp_ave)
    
    scale_rows <- function(x) t(scale(t(x)))
    exp_ave_scale <- scale_rows(exp_ave)
>>>>>>> f40c6eb162999aea48b4d3f321b5e2bce1953d15
    exp_ave_scale[exp_ave_scale > 2] <- 2
    exp_ave_scale[exp_ave_scale < -2] <- -2
    
    output$heatmapPlot <- renderImage({
<<<<<<< HEAD
      outfile <- tempfile(fileext = ".png")
      png(outfile, width = 1200, height = 900, res = 150)
=======
      # Create a temporary file to save the heatmap
      outfile <- tempfile(fileext = ".png")
      
      png(outfile, width = 1200, height = 900, res = 150)  # puedes ajustar tamaño y resolución
      
      # Dibujar el heatmap (usa tus datos ya procesados)
>>>>>>> f40c6eb162999aea48b4d3f321b5e2bce1953d15
      pheatmap::pheatmap(exp_ave_scale,
                         cluster_rows = FALSE,
                         cluster_cols = TRUE,
                         show_rownames = TRUE,
                         show_colnames = TRUE,
                         fontsize_col = 10,
                         angle_col = 45,
                         treeheight_row = 0,
                         treeheight_col = 0,
                         fontsize = 12,
                         border_color = "white")
<<<<<<< HEAD
      dev.off()
      list(src = outfile, contentType = "image/png", width = "100%", height = "auto", alt = "Heatmap")
=======
      
      dev.off()
      
      list(
        src = outfile,
        contentType = "image/png",
        width = "100%",   # o píxeles si prefieres
        height = "auto",
        alt = "Heatmap"
      )
>>>>>>> f40c6eb162999aea48b4d3f321b5e2bce1953d15
    }, deleteFile = TRUE)
  })
  
  #server UMAPs#####

  
  
  output$umapPlot <- renderPlot({
    
    # Variable seleccionada para colorear
    var_color <- input$umap_color_by
    
    color_pal <- c(
      "0" = "#CDB4DB",
      "1" = "#F78181",
      "2" = "#82CFFD",
      "6" = "#FFB347",
      "4" = "#AEB6BF",
      "5" = "#F7ACB4",
      "7" = "#7FFFD4",
      "10" = "#FFA07A",
      "8" = "#F4D03F",
      "9" = "#CDB4DB",
      "3" = "#8FCFD1",
      "11" = "#2ECC71",
      "CTRL" = "#9932CC",
      "SSc" = "#228B22",
      "Tn"= "#a9dfbb",
      "Tcm"= "#6cd0c4",
      "Treg"= "#fc6daf",
      "Tem"= "#3faccd",
      "CTL"= "#ff0000",
      "Th17"= '#FF55A3',
      "Th2L"= '#FF9F55',
      "Thf"= '#559FFF',
      "Th1_GZMK"= '#55FF9F',
      "Th1_GZMB"= '#FFC255'
    )
    
    meta$clusters <- as.character(meta$R09)
    
    umap_data <- meta %>%
      mutate(var_to_plot = .[[var_color]])  # crea columna auxiliar para aes
    
    # Determinar si es variable continua
    is_cont <- is.numeric(umap_data$var_to_plot)
    
    p <- ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, color = var_to_plot)) +
      geom_point(size = 2, alpha = 0.8) +
      labs(title = paste("UMAP colored by", var_color), color = "") +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.text = element_blank(),
        axis.title = element_text(size = 12),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", size = 1),
        legend.position = "right",
        legend.text = element_text(size = 10)
      )
    
    vars_with_custom_palette <- c("clusters", "condition", "cell_type","cell_type_Th")
    
    if (is_cont) {
      
      p <- p + scale_color_viridis_b(
        option = "plasma",
        direction = -1,
        na.value = "grey80"
      )
      
    } else {
      
      if (var_color %in% vars_with_custom_palette) {
        
        # Usar paleta definida por el usuario
        p <- p + scale_color_manual(
          values = color_pal,
          na.value = "grey80"
        )
        
      } else {
        
        # Paleta discreta automática (fallback)
        levs <- sort(unique(umap_data$var_to_plot))
        n_colors <- length(levs)
        
        palette <- suppressWarnings(
          RColorBrewer::brewer.pal(min(max(3, n_colors), 8), "Set2")
        )
        
        if (n_colors > 8) {
          palette <- colorRampPalette(palette)(n_colors)
        }
        
        p <- p + scale_color_manual(
          values = setNames(palette, levs),
          na.value = "grey80"
        )
      }
    }
    
    p <- p + coord_fixed() 
    p
  })

  
  #server UMAPs gene expr#####
<<<<<<< HEAD

=======
  
>>>>>>> f40c6eb162999aea48b4d3f321b5e2bce1953d15
  observe({
    updateSelectizeInput(
      session,
      inputId = "gene_umap",
      choices = rownames(exprs_norm),
      server = TRUE
    )
  })
  
  output$umapPlot_expression <- renderPlot({
    req(input$gene_umap)
    FeaPlot(
      geneName = input$gene_umap,
      meta = meta,
      expres_norm = exprs_norm
    )
  })
  
  
  #server Boxplot expr#####
  
<<<<<<< HEAD
  observe({
    updateSelectizeInput(session, "gene_boxplot",
                         choices = rownames(exprs_norm),
                         selected = "CD4",
                         server = TRUE)
  })
  
  # Individual-level boxplot
  output$boxplot_expr_individual <- renderPlot({
    req(input$gene_boxplot, input$group_by)
    plot_boxplot_expr_individual(
      gene = input$gene_boxplot,
      group_var = input$group_by,
      exprs_norm = exprs_norm,
      meta = meta
    )
  })
  
  # Grouped boxplot
  output$boxplot_expr_grouped <- renderPlot({
    req(input$gene_boxplot, input$group_by)
    plot_boxplot_expr_grouped(
      gene = input$gene_boxplot,
      group_var = input$group_by,
      exprs_norm = exprs_norm,
      meta = meta
    )
=======
  #     meta$clusters <- as.character(meta$R09)
  
  # 1. Per-sample individual expression (como al principio)
  output$boxplot_expr_individual <- renderPlot({
    req(input$gene_boxplot, input$group_by)
    
    group_var <- input$group_by
    
    # Extraer expresión del gen
    gene_expr <- exprs_norm[input$gene_boxplot, ]
    
    meta$clusters <- as.character(meta$R09)
    
    # Juntar con metadata
    df <- meta %>%
      mutate(expression = gene_expr[colnames(exprs_norm)]) %>%
      group_by(Sample, .data[[group_var]]) %>%
      summarise(mean_expr = mean(expression, na.rm = TRUE), .groups = "drop")
    
    # Definir colores (puedes ajustar esta parte para más variables)
    custom_colors <- c(
      "CTRL" = "#D8BFD8", "SSc" = "#C3E5C3",
      "Tn"= "#a9dfbb", "Tcm"= "#6cd0c4", "Treg"= "#fc6daf",
      "Tem"= "#3faccd", "CTL"= "#ff0000",
      "0" = "#CDB4DB", "1" = "#F78181", "2" = "#82CFFD", "3" = "#8FCFD1",
      "4" = "#AEB6BF", "5" = "#F7ACB4", "6" = "#FFB347", "7" = "#7FFFD4",
      "8" = "#F4D03F", "9" = "#CDB4DB", "10" = "#FFA07A", "11" = "#2ECC71"
    )
    
    # Forzar orden de niveles según la variable agrupadora
    if (group_var == "cell_type") {
      df[[group_var]] <- factor(df[[group_var]], levels = c("Tn", "Tcm", "Tem", "CTL", "Treg"))
    } else if (group_var == "clusters") {
      df[[group_var]] <- factor(df[[group_var]], levels = as.character(0:11))
    } else if (group_var == "condition") {
      df[[group_var]] <- factor(df[[group_var]], levels = c("CTRL", "SSc"))
    }
    
    p <- ggplot(df, aes(x = .data[[group_var]],
                        y = mean_expr,
                        fill = .data[[group_var]])) +
      geom_boxplot(width = 0.5, alpha = 0.75, outlier.shape = NA) +
      geom_point(alpha = 0.9,
                 shape = 21,
                 colour = "black",
                 stroke = 0.3,
                 size = 2.5,
                 position = position_jitter(width = 0.15)) +
      labs(title = paste("Expression of", input$gene_boxplot),
           x = "",
           y = "Mean expression per sample",
           fill = "") +
      scale_fill_manual(values = custom_colors, na.translate = FALSE) +
      theme_bw() +
      theme(legend.position = 'none',
            axis.text.x = element_text(size = 12, face = "bold"),
            axis.text.y = element_text(size = 12),
            axis.title.y = element_text(size = 20),
            title = element_text(size = 18),
            plot.margin = unit(c(1, 1, 1, 1), "mm"))
    
    # Comparación estadística si solo hay dos grupos
    #if (requireNamespace("ggpubr", quietly = TRUE) && length(unique(df[[group_var]])) == 2) {
      #p <- p + ggpubr::stat_compare_means(
        #comparisons = list(unique(df[[group_var]])),
        #label = "p.signif",
        #method = "wilcox.test"
      #)
    #}
    
    p
  })
  
  # 2. Agrupado por variable (e.g. R09 o cell_type)
  output$boxplot_expr_grouped <- renderPlot({
    req(input$gene_boxplot, input$group_by)
    
    group_var <- input$group_by
    gene_expr <- exprs_norm[input$gene_boxplot, ]
    
    meta$clusters <- as.character(meta$R09)
    
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
      labs(title = paste("Expression of", input$gene_boxplot, "(grouped by", group_var, ")"),
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
>>>>>>> f40c6eb162999aea48b4d3f321b5e2bce1953d15
  })
  
  # server GSEA ####
  gsea_names_filtered <- reactive({
    req(input$gsea_category)
    
    suffix <- switch(input$gsea_category,
                     "H"  = "_Hallmark",
<<<<<<< HEAD
                     "Reactome" = "_C2")
=======
                     "C2" = "_C2")
>>>>>>> f40c6eb162999aea48b4d3f321b5e2bce1953d15
    
    names(gseaTF_list)[endsWith(names(gseaTF_list), suffix)]
  })
  
  observeEvent(input$gsea_category, {
<<<<<<< HEAD
    choices <- gsea_names_filtered()
    req(length(choices) > 0)
    
    selected_value <- if (input$gsea_category == "H" && "Global_Hallmark" %in% choices) {
      "Global_Hallmark"
    } else {
      choices["Global_Reactome"]
    }
=======
    
    choices <- gsea_names_filtered()
    
    req(length(choices) > 0)
>>>>>>> f40c6eb162999aea48b4d3f321b5e2bce1953d15
    
    updateSelectInput(
      session,
      "gsea_dataset",
      choices = choices,
<<<<<<< HEAD
      selected = selected_value
=======
      selected = choices[1]
>>>>>>> f40c6eb162999aea48b4d3f321b5e2bce1953d15
    )
  })
  
  gsea_filtered <- reactive({
    req(input$gsea_dataset)
    df <- gseaTF_list[[input$gsea_dataset]]
    req(df)
    
    df %>%
      filter(padj < input$fdr_gsea) %>%
      arrange(desc(NES))
  })
  
  output$gsea_table <- renderDT({
    datatable(gsea_filtered(), options = list(pageLength = 20), rownames = FALSE)
  })
  
  output$gsea_plot <- renderPlot({
    df <- gsea_filtered()
    
    # Filtrar top up y down
    top_pos <- df %>% filter(NES > 0) %>% top_n(input$top_n_gsea, NES)
    top_neg <- df %>% filter(NES < 0) %>% top_n(-input$top_n_gsea, NES)
    top_df <- bind_rows(top_pos, top_neg)
    
    ggplot(top_df, aes(x = reorder(pathway, NES), y = NES, fill = NES > 0)) +
      geom_col(aes(fill = ifelse(NES > 0, "NES > 0", "NES < 0")), alpha = 0.8,
               width = 0.75, colour = "black") +
      coord_flip() +
      scale_fill_manual(values = c("NES > 0" = "#F77663", "NES < 0" = "#3D799E")) +
      labs(x = "", y = "Normalized Enrichment Score (NES)", fill = "Direction") +
      theme_bw() +
      theme(
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 14),
        legend.position = "top",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
  })

  
  
}

# Run the application 
shinyApp(ui = ui, server = server)

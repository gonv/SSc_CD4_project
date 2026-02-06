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


# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # theme
  theme = bslib::bs_theme(version = 4, bootswatch = "litera"),
  
  # Logos
  column(
    12,
    tags$div(
      style = "display: flex; justify-content: center; align-items: center; gap: 40px; padding: 10px;",
      tags$img(src = "MainPage/project_logo.png", style = "height:170px;"),
      tags$img(src = "MainPage/ipblnlogo.jpg", style = "height:60px;"),
      tags$img(src = "MainPage/BCH-logo.png", style = "height:70px;"),
      tags$img(src = "MainPage/Broad-logo.jpeg", style = "height:40px;"),
      tags$img(src = "MainPage/MGA_logo_v6.png", style = "height:100px;")
    )
  ),

    # Application title
    #titlePanel("SSc CD4 project"),
  
  
 tabsetPanel(
    
   # ui Home ####
   
   tabPanel("Home",
            fluidRow(
              column(12,
                     h2("Welcome to the CD4+ T cells in SSc resource"),
                     br(),
                     
                     h4("About this app"),
                     p("This Shiny App provides an interactive exploration of transcriptional signatures in CD4+ T cells from patients with Systemic Sclerosis (SSc). 
                    It includes differential expression analyses, gene set enrichment results, clustering information, and cell-type specific markers."),
                     
                     p("For more details, please refer to the preprint:"),
                     tags$a(href = "https://www.biorxiv.org/content/10.64898/2026.01.03.697349v1", 
                            "Preprint on bioRxiv", target = "_blank"),
                     
                     br(),
                     
                     p("The full code and datasets used to build this application are available at:"),
                     tags$a(href = "https://github.com/gonv/SSc_CD4_project", 
                            "GitHub Repository", target = "_blank"),
                     
                     br(), br(),
                     
                     h4("Citation"),
                     p("If you use this resource, please cite:"),
                     tags$code('Villanueva-Martín G, et al. (2026). Transcriptional landscape of CD4+ T cells in Systemic Sclerosis. bioRxiv.')
              )
            )
   ),
   
   # ui demographic ####
   tabPanel("Demographics",
            fluidRow(
              column(12,
                     h2("Patient Demographics"),
                     DT::DTOutput("demographic_table")
              )
            )
   ),
    
    #ui Table####
    tabPanel("Differential Expression analysis",

             sidebarLayout(
               sidebarPanel(
                 fluidRow(
                   column(12,
                          h4("Select Differential expression analysis SSc vs CTRL"),
                          selectInput("dea_source", "",
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
                        h4("Volcano plot (Global DEA only)")
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
                 sliderInput("pval_cutoff_volcano", "FDR threshold:", min = 0, max = 1, value = 0.05, step = 0.01),
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
                                    selected = "Global",
                                    width = "auto"),
                        selectizeInput("genes_heatmap_only", 
                                       "Genes to plot:",
                                       choices = NULL,
                                       multiple = TRUE,
                                       selected = NULL,
                                       options = list(maxItems = 50),
                                       width = "auto"),
                        actionButton("runHeatmap", "Run Heatmap",
                                     style = "width: 10ß0%; padding: 6px; font-size: 16px;")
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
               sidebarPanel(width = 12,
                 selectInput("umap_color_by", "Color by:",
                             choices = c("clusters", "cell_type", "condition", "Sample",
                                         "percent_mito"), #"cellType_Th",
                             selected = "clusters")
               ),
               mainPanel(width = 12,
                 plotOutput("umapPlot", height = "600px")
               )
             )
    ),
    
    
    #ui UMAP expression#####
    
    tabPanel("Gene expression in UMAP",
             sidebarLayout(
               sidebarPanel(width = 12,
                 selectizeInput(
                   options = list(placeholder = 'Type a gene name...'),
                   inputId = "gene_umap",
                   label = "Select gene:",
                   selected = "CD4",
                   choices = NULL,
                   multiple = FALSE
                 )
               ),
               mainPanel(width = 12,
                 plotOutput("umapPlot_expression", height = "600px")
               )
             )
    ),
    

    #ui Box plot#####
    
    tabPanel("Gene Expression Boxplot",
             fluidRow(
               column(6,
                      selectizeInput("gene_boxplot", "Select gene:", 
                                     choices = NULL,  
                                     selected = "CD4",
                                     options = list(maxOptions = 5000))
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
                                   choices = c("Hallmark" = "H", "Reactome" = "C2"),
                                   selected = "H", inline = TRUE)
               ),
               column(4,
                      selectInput("gsea_dataset", "Select GSEA dataset:",
                                  choices = names(gseaTF_list), 
                                  selected = if ("Global_Hallmark" %in% names(gseaTF_list)) "Global_Hallmark" else "Global_Reactome" %in% names(gseaTF_list))
               ),
               column(4,
                      numericInput("fdr_gsea", "FDR cutoff:", value = 0.1, step = 0.01),
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
  

  demographic_info <- read.delim("data/demographic_info.txt", sep = "\t", header = TRUE)
  
  # server demographic ####
  
  output$demographic_table <- DT::renderDT({
    DT::datatable(demographic_info, 
                  options = list(pageLength = 10, autoWidth = TRUE),
                  rownames = FALSE)
  })
  
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
                filter = 'top')
    })
  })
  
  
  #server Volcano #####
  
  observe({
    updateSelectizeInput(
      session,
      inputId = "genes_heatmap_volcano",
      choices = unique(dea_list[["Global"]]$gene),
      server = TRUE
    )
  })
  
  observeEvent(input$runVolcanoOnly, {
    
    volcano_data <- dea_list[["Global"]]
    
    volcano_data$significance <- case_when(
      volcano_data$adj.P.Val < input$pval_cutoff_volcano & volcano_data$Estimate > 0 ~ "Upregulated",
      volcano_data$adj.P.Val < input$pval_cutoff_volcano & volcano_data$Estimate < 0 ~ "Downregulated",
      TRUE ~ "Not_significant"
    )
    
    volcano_data$delabel <- ifelse(
      volcano_data$gene %in% input$genes_heatmap_volcano,
      volcano_data$gene,
      NA
    )
    
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
        geom_hline(yintercept = fdr_threshold, linetype = "dashed", color = "blue") +
        geom_vline(xintercept = 0, linetype = "dashed", color = "black")
    })
  })
  
  #server heatmap #####
  
  #exprs_matrix <- reactive({
    #readRDS("data/exprs_norm.rds")
  #})
  
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
    req(input$dea_source_heatmap)
    
    dea_data <- dea_list[[input$dea_source_heatmap]]
    
    if (length(input$genes_heatmap_only) > 0) {
      top_genes <- input$genes_heatmap_only
    } else {
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
    exp_ave_scale[exp_ave_scale > 2] <- 2
    exp_ave_scale[exp_ave_scale < -2] <- -2
    
    output$heatmapPlot <- renderImage({
      outfile <- tempfile(fileext = ".png")
      png(outfile, width = 1200, height = 900, res = 150)
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
      dev.off()
      list(src = outfile, contentType = "image/png", width = "100%", height = "auto", alt = "Heatmap")
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
  })
  
  # server GSEA ####
  gsea_names_filtered <- reactive({
    req(input$gsea_category)
    
    suffix <- switch(input$gsea_category,
                     "H"  = "_Hallmark",
                     "C2" = "_Reactome")
    
    names(gseaTF_list)[endsWith(names(gseaTF_list), suffix)]
  })
  
  observeEvent(input$gsea_category, {
    choices <- gsea_names_filtered()
    req(length(choices) > 0)
    
    selected_value <- if (input$gsea_category == "H" && "Global_Hallmark" %in% choices) {
      "Global_Hallmark"
    } else {
      choices["Global_Reactome"]
    }
    
    updateSelectInput(
      session,
      "gsea_dataset",
      choices = choices,
      selected = selected_value
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

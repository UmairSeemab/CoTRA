# =======================================================================
# modules/sc/variable_pca.R
# CoTRA scRNA-seq: Variable Genes, HVGs, PCA visualization
# With unified color palette + final_annotation
# =======================================================================

library(shiny)
library(Seurat)
library(dplyr)
library(ggplot2)
library(plotly)
library(corrplot)

utils::globalVariables("GvHD")

mod_sc_variable_pca_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    fluidRow(
      column(
        6,
        conditionalPanel(
          condition = "output.check_vg",
          ns = ns,
          bs4Dash::bs4Card(
            title = tagList("PCA Parameters",
                            actionButton(ns("HELP_PCA_PARAMS"), "", icon = icon("question-circle"), class = "btn-xs")),
            width = 12,
            solidHeader = TRUE,
            status = "primary",
            collapsible = TRUE,
            uiOutput(ns("list_vg")),
            uiOutput(ns("nb_pca")),
            textInput(ns("pca_name"), "Enter a Name for the PCA Result", "pca", width = "50%"),
            uiOutput(ns("select_regres")),
            actionButton(ns("run_pca"), label = "Run PCA")
          )
        ),
        conditionalPanel(
          condition = "output.check_pca",
          ns = ns,
          bs4Dash::bs4Card(
            title = tagList("PC versus Metadata Scatter",
                            actionButton(ns("HELP_PCA_SCATTER"), "", icon = icon("question-circle"), class = "btn-xs")),
            width = 12,
            solidHeader = TRUE,
            status = "primary",
            collapsible = TRUE,
            plotlyOutput(ns("plotPCA")),
            fluidRow(
              column(6, uiOutput(ns("select_PCA"))),
              column(6, uiOutput(ns("select_numeric")))
            ),
            hr(),
            radioButtons(
              ns("pca_scatter_format"),
              "Download format",
              choices = c("PNG" = "png", "SVG" = "svg"),
              inline = TRUE
            ),
            downloadButton(ns("down_pca_scatter"), "Download", class = "btn btn-primary")
          )
        )
      ),
      column(
        6,
        conditionalPanel(
          condition = "output.check_pca",
          ns = ns,
          bs4Dash::bs4Card(
            title = tagList("Scree Plot",
                            actionButton(ns("HELP_PCA_PC"), "", icon = icon("question-circle"), class = "btn-xs")),
            width = 12,
            solidHeader = TRUE,
            status = "primary",
            collapsible = TRUE,
            plotlyOutput(ns("ScreePlot"), height = 300),
            hr(),
            radioButtons(
              ns("scree_format"),
              "Download format",
              choices = c("PNG" = "png", "SVG" = "svg"),
              inline = TRUE
            ),
            downloadButton(ns("down_scree"), "Download", class = "btn btn-primary")
          ),
          bs4Dash::bs4Card(
            title = tagList("PC versus Metadata Correlation Heatmap",
                            actionButton(ns("HELP_PCA_HEATMAP"), "", icon = icon("question-circle"), class = "btn-xs")),
            width = 12,
            solidHeader = TRUE,
            status = "primary",
            collapsible = TRUE,
            plotlyOutput(ns("CorHeat"), height = 350),
            hr(),
            radioButtons(
              ns("corheat_format"),
              "Download format",
              choices = c("PNG" = "png", "SVG" = "svg"),
              inline = TRUE
            ),
            downloadButton(ns("down_corheat"), "Download", class = "btn btn-primary")
          ),
          bs4Dash::bs4Card(
            title = "PC Loadings Heatmap (DimHeatmap)",
            width = 12,
            solidHeader = TRUE,
            status = "primary",
            collapsible = TRUE,
            uiOutput(ns("select_PC_heat")),
            plotOutput(ns("DimHeatmap"), height = 300),
            hr(),
            radioButtons(
              ns("dimheat_format"),
              "Download format",
              choices = c("PNG" = "png", "SVG" = "svg"),
              inline = TRUE
            ),
            downloadButton(ns("down_DimHeatmap"), "Download", class = "btn btn-primary")
          )
        )
      )
    )
  )
}

mod_sc_variable_pca_server <- function(id, rval) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    params <- reactiveValues(
      pca_name = NULL,
      current_seu_name = NULL
    )
    
    observeEvent(input[["HELP_PCA_PARAMS"]], {
      showModal(modalDialog(
        title = HTML(help_infos[help_infos$key == "HELP_PCA_PARAMS", "title"]),
        HTML(help_infos[help_infos$key == "HELP_PCA_PARAMS", "value"]),
        easyClose = TRUE,
        footer = modalButton("OK")
      ))
    })
    
    observeEvent(input[["HELP_PCA_SCATTER"]], {
      showModal(modalDialog(
        title = HTML(help_infos[help_infos$key == "HELP_PCA_SCATTER", "title"]),
        HTML(help_infos[help_infos$key == "HELP_PCA_SCATTER", "value"]),
        easyClose = TRUE,
        footer = modalButton("OK")
      ))
    })
    
    observeEvent(input[["HELP_PCA_PC"]], {
      showModal(modalDialog(
        title = HTML(help_infos[help_infos$key == "HELP_PCA_PC", "title"]),
        HTML(help_infos[help_infos$key == "HELP_PCA_PC", "value"]),
        easyClose = TRUE,
        footer = modalButton("OK")
      ))
    })
    
    observeEvent(input[["HELP_PCA_HEATMAP"]], {
      showModal(modalDialog(
        title = HTML(help_infos[help_infos$key == "HELP_PCA_HEATMAP", "title"]),
        HTML(help_infos[help_infos$key == "HELP_PCA_HEATMAP", "value"]),
        easyClose = TRUE,
        footer = modalButton("OK")
      ))
    })
    
    output$list_vg <- renderUI({
      if (is.null(rval$seurat) || is.null(rval$seurat_selected)) return(NULL)
      ns <- session$ns
      
      gene_lists <- character(0)
      if (!is.null(rval$genes_list[[rval$seurat_selected]])) {
        gene_lists <- c(gene_lists, names(rval$genes_list[[rval$seurat_selected]]))
      }
      if (!is.null(rval$genes_list[["imported"]])) {
        gene_lists <- c(gene_lists, names(rval$genes_list[["imported"]]))
      }
      gene_lists <- setdiff(gene_lists, c("s.genes", "g2m.genes"))
      
      if (length(gene_lists) == 0) {
        params$pca_name <- NULL
        return(NULL)
      }
      
      selectInput(
        ns("list_vg"),
        label = "Choose Gene List",
        choices = gene_lists,
        selected = tail(gene_lists, 1)
      )
    })
    
    output$nb_pca <- renderUI({
      if (is.null(rval$seurat)) return(NULL)
      ns <- session$ns
      numericInput(
        ns("nb_pca"),
        label = "Set the Number of PCs to Compute",
        value = 50,
        min = 10,
        max = 1000,
        step = 1,
        width = "50%"
      )
    })
    
    numericCol <- reactive({
      if (is.null(rval$seurat)) return(NULL)
      colnames(dplyr::select_if(rval$seurat@meta.data, is.numeric))
    })
    
    output$select_regres <- renderUI({
      if (is.null(rval$seurat)) return(NULL)
      ns <- session$ns
      nums <- numericCol()
      if (is.null(nums) || length(nums) == 0) return(NULL)
      selectizeInput(
        ns("to_regres"),
        label = "Select Variable(s) to Regress Out",
        choices = nums,
        selected = NULL,
        multiple = TRUE,
        width = "50%"
      )
    })
    
    observeEvent(input$run_pca, {
      if (is.null(rval$seurat) || is.null(rval$seurat_selected)) return(NULL)
      
      if (is.na(input$nb_pca) || input$nb_pca < 10 || input$nb_pca > 1000) {
        showModal(modalDialog(
          "The number of PCs must be between 10 and 1000.",
          easyClose = TRUE,
          footer = modalButton("OK")
        ))
        return(NULL)
      }
      
      if (input$pca_name == "") {
        showModal(modalDialog(
          "Please enter a PCA name.",
          easyClose = TRUE,
          footer = modalButton("OK")
        ))
        return(NULL)
      }
      
      if (!grepl("^pca", input$pca_name, ignore.case = TRUE) ||
          grepl("[^A-Za-z0-9_]", input$pca_name)) {
        showModal(modalDialog(
          "PCA name must start with 'pca' and contain only letters, numbers, or underscore.",
          easyClose = TRUE,
          footer = modalButton("OK")
        ))
        return(NULL)
      }
      
      if (input$pca_name %in% Reductions(rval$seurat)) {
        showModal(modalDialog(
          "This PCA name already exists.",
          easyClose = TRUE,
          footer = modalButton("OK")
        ))
        return(NULL)
      }
      
      if (is.null(input$list_vg)) {
        showModal(modalDialog(
          "No gene list selected.",
          easyClose = TRUE,
          footer = modalButton("OK")
        ))
        return(NULL)
      }
      
      show_modal_spinner(text = "Running regression and PCA, please wait...", spin = "circle")
      
      all.genes <- rownames(rval$seurat)
      vars_to_regress <- input$to_regres
      
      rval$seurat <- ScaleData(
        rval$seurat,
        features = all.genes,
        do.scale = FALSE,
        do.center = TRUE,
        vars.to.regress = vars_to_regress
      )
      
      if (input$list_vg %in% names(rval$genes_list[[rval$seurat_selected]])) {
        gene_list <- rval$genes_list[[rval$seurat_selected]][[input$list_vg]]
      } else if (!is.null(rval$genes_list[["imported"]]) &&
                 input$list_vg %in% names(rval$genes_list[["imported"]])) {
        gene_list <- rval$genes_list[["imported"]][[input$list_vg]]
      } else {
        gene_list <- NULL
      }
      
      if (is.null(gene_list) || length(gene_list) < 10) {
        remove_modal_spinner()
        showModal(modalDialog(
          "Selected gene list is empty or too small for PCA.",
          easyClose = TRUE,
          footer = modalButton("OK")
        ))
        return(NULL)
      }
      
      rval$seurat <- RunPCA(
        rval$seurat,
        features = gene_list,
        npcs = input$nb_pca,
        reduction.name = input$pca_name,
        assay = if (!is.null(rval$norm_assay)) rval$norm_assay else DefaultAssay(rval$seurat)
      )
      
      if (is.null(rval$parameters[[rval$seurat_selected]])) {
        rval$parameters[[rval$seurat_selected]] <- list()
      }
      
      rval$parameters[[rval$seurat_selected]][[input$pca_name]] <- c(
        genes = input$list_vg,
        npcs = input$nb_pca,
        regress = if (is.null(vars_to_regress)) NA else paste(vars_to_regress, collapse = ",")
      )
      
      if (is.null(rval$parameters[[rval$seurat_selected]][["steps"]])) {
        rval$parameters[[rval$seurat_selected]][["steps"]] <- "pca"
      } else {
        rval$parameters[[rval$seurat_selected]][["steps"]] <- c(
          rval$parameters[[rval$seurat_selected]][["steps"]],
          "pca"
        )
      }
      
      params$pca_name <- input$pca_name
      params$current_seu_name <- rval$seurat_selected
      
      if (!is.null(rval$output)) {
        saveRDS(
          reactiveValuesToList(rval),
          file.path(rval$output, "Project_elements.rds")
        )
      }
      
      remove_modal_spinner()
    })
    
    nbPCACol <- reactive({
      if (is.null(params$pca_name) || is.null(rval$seurat)) return(NULL)
      colnames(Embeddings(rval$seurat, reduction = params$pca_name))
    })
    
    output$select_PCA <- renderUI({
      if (is.null(params$pca_name) || is.null(rval$seurat)) return(NULL)
      ns <- session$ns
      pcs <- nbPCACol()
      if (is.null(pcs)) return(NULL)
      selectInput(
        ns("axisPCA"),
        "X Axis (PC)",
        choices = pcs,
        selected = pcs[1]
      )
    })
    
    output$select_numeric <- renderUI({
      if (is.null(params$pca_name) || is.null(rval$seurat)) return(NULL)
      ns <- session$ns
      nums <- numericCol()
      if (is.null(nums) || length(nums) == 0) return(NULL)
      selectInput(
        ns("numeric"),
        "Y Axis (numeric metadata)",
        choices = nums,
        selected = nums[1]
      )
    })
    
    pca_scatter_data <- reactive({
      if (is.null(rval$seurat) || is.null(input$axisPCA) || is.null(input$numeric)) return(NULL)
      vars <- c(input$axisPCA, input$numeric)
      vars <- unique(c(vars, "Project"))
      df <- FetchData(rval$seurat, vars = vars)
      df
    })
    
    pca_scatter_gg <- reactive({
      df <- pca_scatter_data()
      if (is.null(df)) return(NULL)
      x_var <- colnames(df)[1]
      y_var <- colnames(df)[2]
      cor_val <- suppressWarnings(cor(df[[x_var]], df[[y_var]], method = "spearman", use = "complete.obs"))
      g <- ggplot(df, aes_string(x = x_var, y = y_var, color = "Project")) +
        geom_point(alpha = 0.7, size = 0.6) +
        theme_minimal() +
        theme(legend.position = "right") +
        labs(
          x = input$axisPCA,
          y = input$numeric,
          title = paste0("Spearman correlation: ", round(cor_val, 2))
        )
      g
    })
    
    output$plotPCA <- renderPlotly({
      g <- pca_scatter_gg()
      if (is.null(g)) return(NULL)
      ggplotly(g)
    })
    
    output$down_pca_scatter <- downloadHandler(
      filename = function() {
        paste0("pca_metadata_scatter.", input$pca_scatter_format)
      },
      content = function(file) {
        g <- pca_scatter_gg()
        if (is.null(g)) return(NULL)
        if (input$pca_scatter_format == "png") {
          png(file, width = 1800, height = 1500, res = 300)
          print(g)
          dev.off()
        } else {
          svg(file, width = 7, height = 6)
          print(g)
          dev.off()
        }
      }
    )
    
    scree_df <- reactive({
      if (is.null(params$pca_name) || is.null(rval$seurat)) return(NULL)
      emb <- Embeddings(rval$seurat, reduction = params$pca_name)
      if (is.null(emb) || ncol(emb) == 0) return(NULL)
      vars <- apply(emb, 2, var)
      pcs <- seq_along(vars)
      df <- data.frame(
        PC = pcs,
        Variance = vars,
        VarExplained = vars / sum(vars)
      )
      df
    })
    
    scree_gg <- reactive({
      df <- scree_df()
      if (is.null(df)) return(NULL)
      g <- ggplot(df, aes(x = PC, y = VarExplained)) +
        geom_line() +
        geom_point() +
        theme_minimal() +
        labs(
          x = "Principal Component",
          y = "Proportion of variance explained"
        )
      g
    })
    
    output$ScreePlot <- renderPlotly({
      g <- scree_gg()
      if (is.null(g)) return(NULL)
      ggplotly(g)
    })
    
    output$down_scree <- downloadHandler(
      filename = function() {
        paste0("pca_scree.", input$scree_format)
      },
      content = function(file) {
        g <- scree_gg()
        if (is.null(g)) return(NULL)
        if (input$scree_format == "png") {
          png(file, width = 1800, height = 1500, res = 300)
          print(g)
          dev.off()
        } else {
          svg(file, width = 7, height = 6)
          print(g)
          dev.off()
        }
      }
    )
    
    dfHeatCor <- reactive({
      if (is.null(params$pca_name) || is.null(rval$seurat)) return(NULL)
      pcs <- nbPCACol()
      if (is.null(pcs)) return(NULL)
      ndims <- min(length(pcs), 20)
      emb <- Embeddings(rval$seurat, reduction = params$pca_name)[, 1:ndims, drop = FALSE]
      nums <- numericCol()
      if (is.null(nums) || length(nums) == 0) return(NULL)
      df <- cbind(emb, rval$seurat@meta.data[, nums, drop = FALSE])
      cmat <- suppressWarnings(cor(df, method = "spearman", use = "pairwise.complete.obs"))
      cmat <- cmat[nums, colnames(emb), drop = FALSE]
      cmat
    })
    
    corheat_plotly <- reactive({
      cmat <- dfHeatCor()
      if (is.null(cmat)) return(NULL)
      plot_ly(
        x = colnames(cmat),
        y = rownames(cmat),
        z = cmat,
        type = "heatmap",
        colors = grDevices::colorRamp(c("#313695", "#ABD9E9", "#FFFFBF", "#FDAE61", "#A50026"))
      ) %>%
        layout(
          xaxis = list(title = "PCs"),
          yaxis = list(title = "Metadata", autorange = "reversed")
        )
    })
    
    output$CorHeat <- renderPlotly({
      corheat_plotly()
    })
    
    output$down_corheat <- downloadHandler(
      filename = function() {
        paste0("pca_metadata_correlation.", input$corheat_format)
      },
      content = function(file) {
        cmat <- dfHeatCor()
        if (is.null(cmat)) return(NULL)
        if (input$corheat_format == "png") {
          png(file, width = 1800, height = 1800, res = 300)
          corrplot::corrplot(cmat, tl.col = "black", tl.srt = 45, type = "full")
          dev.off()
        } else {
          svg(file, width = 7, height = 7)
          corrplot::corrplot(cmat, tl.col = "black", tl.srt = 45, type = "full")
          dev.off()
        }
      }
    )
    
    output$select_PC_heat <- renderUI({
      if (is.null(params$pca_name)) return(NULL)
      ns <- session$ns
      pcs <- nbPCACol()
      if (is.null(pcs)) return(NULL)
      idx <- seq_along(pcs)
      selectInput(
        ns("PCheat"),
        "PC for DimHeatmap",
        choices = idx,
        selected = idx[1],
        width = "25%"
      )
    })
    
    DimHeatmap_plot <- reactive({
      if (is.null(params$pca_name) || is.null(rval$seurat) || is.null(input$PCheat)) return(NULL)
      DimHeatmap(
        rval$seurat,
        reduction = params$pca_name,
        dims = as.numeric(input$PCheat),
        cells = 500,
        balanced = TRUE
      )
    })
    
    output$DimHeatmap <- renderPlot({
      DimHeatmap_plot()
    })
    
    output$down_DimHeatmap <- downloadHandler(
      filename = function() {
        paste0("pca_dimheatmap.", input$dimheat_format)
      },
      content = function(file) {
        if (input$dimheat_format == "png") {
          png(file, width = 1800, height = 1800, res = 300)
          if (!is.null(params$pca_name) && !is.null(input$PCheat)) {
            DimHeatmap(
              rval$seurat,
              reduction = params$pca_name,
              dims = as.numeric(input$PCheat),
              cells = 500,
              balanced = TRUE
            )
          }
          dev.off()
        } else {
          svg(file, width = 7, height = 7)
          if (!is.null(params$pca_name) && !is.null(input$PCheat)) {
            DimHeatmap(
              rval$seurat,
              reduction = params$pca_name,
              dims = as.numeric(input$PCheat),
              cells = 500,
              balanced = TRUE
            )
          }
          dev.off()
        }
      }
    )
    
    output$check_pca <- reactive({
      if (is.null(params$current_seu_name)) return(FALSE)
      if (is.null(params$pca_name)) return(FALSE)
      if (is.null(rval$seurat_selected)) return(FALSE)
      params$current_seu_name == rval$seurat_selected && params$pca_name %in% Reductions(rval$seurat)
    })
    outputOptions(output, "check_pca", suspendWhenHidden = FALSE)
    
    output$check_vg <- reactive({
      if (is.null(rval$seurat) || is.null(rval$seurat_selected)) return(FALSE)
      !is.null(rval$genes_list[[rval$seurat_selected]])
    })
    outputOptions(output, "check_vg", suspendWhenHidden = FALSE)
    
    rval
    
    
  })
}
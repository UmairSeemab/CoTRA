# ==========================================================
# modules/sc/demultiplexing.R
# CoTRA scRNA-seq: HTO Demultiplexing module
# ==========================================================

mod_sc_demux_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    fluidRow(
      bs4Dash::bs4Card(
        title = "HTO Counts Normalization",
        width = 3,
        collapsible = TRUE,
        solidHeader = TRUE,
        status = "primary",
        selectInput(
          ns("norm_hto"),
          "Select normalization method",
          c(
            "Centred log ratio" = "clr_m",
            "Log normalization" = "ln_m",
            "Relative counts" = "rc_m"
          )
        ),
        actionButton(ns("go_norm_hto"), label = "Normalize")
      ),
      
      conditionalPanel(
        condition = "output.norm_hto",
        ns = ns,
        bs4Dash::bs4Card(
          title = "Sample demultiplexing using MULTIseqDemux",
          width = 4,
          collapsible = TRUE,
          solidHeader = TRUE,
          status = "primary",
          radioButtons(
            ns("autoThresh_hto"),
            label = "Automatic threshold quantile",
            choices = c("Yes" = "auto_yes", "No" = "auto_no"),
            inline = TRUE
          ),
          uiOutput(ns("quantile_hto")),
          actionButton(ns("MultiseqDemux"), label = "Run demultiplexing")
        )
      )
    ),
    
    fluidRow(
      conditionalPanel(
        condition = "output.is_demux",
        ns = ns,
        
        bs4Dash::bs4Card(
          title = "Demultiplexing results",
          status = "primary",
          solidHeader = TRUE,
          width = 3,
          tableOutput(ns("datatable_demux")),
          actionButton(ns("filter_hto"), label = "Filter out multiplets and negatives")
        ),
        
        conditionalPanel(
          condition = "output.is_UMAP_hto",
          ns = ns,
          bs4Dash::bs4Card(
            title = "HTO data visualization",
            status = "primary",
            solidHeader = TRUE,
            width = 9,
            tabsetPanel(
              tabPanel(
                "UMAP",
                plotOutput(ns("UMAP_hto")),
                br(),
                downPlotUI(id = ns("UMAP_hto_export"))
              ),
              tabPanel(
                "Ridge plots",
                plotOutput(ns("ridge_hto"), height = 650),
                br(),
                downPlotUI(id = ns("ridge_hto_export"))
              )
            )
          )
        )
      )
    )
    
    
  )
}

mod_sc_demux_server <- function(id, seurat) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    rv <- reactiveValues(
      seurat   = NULL,
      quantile = "Auto",
      norm     = NULL
    )
    
    observeEvent(seurat(), {
      if (is.null(rv$seurat)) {
        rv$seurat <- seurat()
      }
    }, ignoreNULL = TRUE)
    
    
    output$quantile_hto <- renderUI({
      if (input$autoThresh_hto == "auto_no") {
        sliderInput(
          ns("quantileman"),
          label = "Quantile to use",
          min   = 0,
          max   = 1,
          value = 0.5,
          step  = 0.01,
          width = "55%"
        )
      } else {
        NULL
      }
    })
    
    
    observeEvent(input$go_norm_hto, {
      req(rv$seurat)
      
      if (!"HTO" %in% Seurat::Assays(rv$seurat)) {
        showModal(
          modalDialog(
            title = "No HTO assay",
            "The current Seurat object does not contain an HTO assay.",
            easyClose = TRUE
          )
        )
        return(NULL)
      }
      
      withProgress(message = "Normalizing HTO counts", value = 0, {
        incProgress(0.2)
        if (input$norm_hto == "clr_m") {
          rv$seurat <- Seurat::NormalizeData(
            rv$seurat,
            normalization.method = "CLR",
            assay       = "HTO",
            scale.factor = 10000,
            margin      = 1,
            verbose     = TRUE
          )
          rv$norm <- "CLR"
        } else if (input$norm_hto == "ln_m") {
          rv$seurat <- Seurat::NormalizeData(
            rv$seurat,
            normalization.method = "LogNormalize",
            assay       = "HTO",
            scale.factor = 10000,
            margin      = 1,
            verbose     = TRUE
          )
          rv$norm <- "Log normalization"
        } else {
          rv$seurat <- Seurat::NormalizeData(
            rv$seurat,
            normalization.method = "RC",
            assay       = "HTO",
            scale.factor = 10000,
            margin      = 1,
            verbose     = TRUE
          )
          rv$norm <- "Relative counts"
        }
        
        incProgress(0.6, detail = "Computing HTO UMAP")
        
        rv$seurat <- Seurat::RunUMAP(
          rv$seurat,
          reduction     = NULL,
          assay         = "HTO",
          features      = rownames(rv$seurat[["HTO"]]),
          reduction.key = "htoUMAP_",
          reduction.name = "hto_UMAP",
          dims          = NULL
        )
        
        incProgress(0.2)
      })
    })
    
    
    observeEvent(input$MultiseqDemux, {
      req(rv$seurat)
      
      if (!"HTO" %in% Seurat::Assays(rv$seurat)) {
        showModal(
          modalDialog(
            title = "No HTO assay",
            "You need an HTO assay to run MULTIseqDemux.",
            easyClose = TRUE
          )
        )
        return(NULL)
      }
      
      withProgress(message = "Running MULTIseqDemux", value = 0, {
        incProgress(0.2)
        
        if (input$autoThresh_hto == "auto_yes") {
          rv$quantile <- "Auto"
          rv$seurat <- Seurat::MULTIseqDemux(
            rv$seurat,
            assay      = "HTO",
            autoThresh = TRUE,
            verbose    = TRUE
          )
        } else {
          if (is.null(input$quantileman) ||
              is.na(input$quantileman) ||
              input$quantileman < 0 ||
              input$quantileman > 1) {
            showModal(
              modalDialog(
                "Quantile must be a number between 0 and 1.",
                easyClose = TRUE
              )
            )
            return(NULL)
          }
          q <- input$quantileman
          rv$quantile <- q
          rv$seurat <- Seurat::MULTIseqDemux(
            rv$seurat,
            assay      = "HTO",
            quantile   = q,
            autoThresh = FALSE,
            verbose    = TRUE
          )
        }
        
        incProgress(0.8)
      })
    })
    
    
    observeEvent(input$filter_hto, {
      req(rv$seurat)
      if (!"MULTI_ID" %in% colnames(rv$seurat@meta.data)) {
        showModal(
          modalDialog(
            title = "No demultiplexing results",
            "Run MULTIseqDemux before filtering.",
            easyClose = TRUE
          )
        )
        return(NULL)
      }
      
      keep_levels <- setdiff(
        as.character(droplevels(rv$seurat@meta.data$MULTI_ID)),
        c("Doublet", "Negative")
      )
      n_cells <- sum(rv$seurat@meta.data$MULTI_ID %in% keep_levels)
      
      if (n_cells < 50) {
        showModal(
          modalDialog(
            title = "Too few cells",
            "Less than 50 cells would remain after filtering. Try a different quantile.",
            easyClose = TRUE
          )
        )
        return(NULL)
      }
      
      showModal(
        modalDialog(
          title = "Filter multiplets and negatives",
          paste0("You will keep ", n_cells, " cells and remove multiplets and negatives."),
          easyClose = FALSE,
          footer = tagList(
            modalButton("Cancel"),
            actionButton(ns("filter_hto_ok"), "OK")
          )
        )
      )
    })
    
    
    observeEvent(input$filter_hto_ok, {
      req(rv$seurat)
      removeModal()
      
      keep_levels <- setdiff(
        as.character(droplevels(rv$seurat@meta.data$MULTI_ID)),
        c("Doublet", "Negative")
      )
      
      Seurat::Idents(rv$seurat) <- "MULTI_ID"
      rv$seurat <- subset(rv$seurat, idents = keep_levels)
    })
    
    
    output$datatable_demux <- renderTable({
      req(rv$seurat)
      if (!"MULTI_ID" %in% colnames(rv$seurat@meta.data)) {
        return(NULL)
      }
      as.data.frame(table(rv$seurat@meta.data[["MULTI_ID"]]))
    })
    
    
    UMAP_hto_plotting <- reactive({
      req(rv$seurat)
      if (!"hto_UMAP" %in% Seurat::Reductions(rv$seurat)) {
        return(NULL)
      }
      if (!"MULTI_ID" %in% colnames(rv$seurat@meta.data)) {
        return(NULL)
      }
      Seurat::DimPlot(
        rv$seurat,
        reduction = "hto_UMAP",
        group.by  = "MULTI_ID"
      ) + ggplot2::ggtitle(paste0("Quantile used: ", rv$quantile))
    })
    
    output$UMAP_hto <- renderPlot({
      UMAP_hto_plotting()
    })
    
    callModule(
      downPlotServer,
      id       = "UMAP_hto_export",
      data     = UMAP_hto_plotting,
      out_file = "UMAP_hto"
    )
    
    
    ridge_hto_plotting <- reactive({
      req(rv$seurat)
      if (!"HTO" %in% Seurat::Assays(rv$seurat)) {
        return(NULL)
      }
      if (!"MULTI_ID" %in% colnames(rv$seurat@meta.data)) {
        return(NULL)
      }
      Seurat::RidgePlot(
        rv$seurat,
        assay   = "HTO",
        features = rownames(rv$seurat[["HTO"]]),
        group.by = "MULTI_ID",
        ncol     = 3
      )
    })
    
    output$ridge_hto <- renderPlot({
      ridge_hto_plotting()
    })
    
    callModule(
      downPlotServer,
      id       = "ridge_hto_export",
      data     = ridge_hto_plotting,
      out_file = "ridge_hto"
    )
    
    
    output$norm_hto <- reactive({
      if (is.null(rv$seurat)) return(FALSE)
      cmds <- Seurat::Command(rv$seurat)
      "NormalizeData.HTO" %in% names(cmds)
    })
    outputOptions(output, "norm_hto", suspendWhenHidden = FALSE)
    
    
    output$is_demux <- reactive({
      if (is.null(rv$seurat)) return(FALSE)
      "MULTI_ID" %in% colnames(rv$seurat@meta.data)
    })
    outputOptions(output, "is_demux", suspendWhenHidden = FALSE)
    
    
    output$is_UMAP_hto <- reactive({
      if (is.null(rv$seurat)) return(FALSE)
      "hto_UMAP" %in% Seurat::Reductions(rv$seurat)
    })
    outputOptions(output, "is_UMAP_hto", suspendWhenHidden = FALSE)
    
    
    list(
      seurat = reactive(rv$seurat)
    )
    
    
  })
}
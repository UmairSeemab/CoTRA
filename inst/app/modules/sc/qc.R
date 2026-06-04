# ==========================================================
# modules/sc/qc.R
# CoTRA scRNA-seq Quality Control Module
#
# Features:
# - Interactive filtering by nFeature_RNA, nCount_RNA,
#   percent.mito, and percent.ribo
# - Slider and numeric threshold control
# - Live retained/removed cell summary
# - Violin and scatter QC plots
# - Downloadable PDF/SVG figures
# - Collapsible help and interpretation sections
# - Safe handling of missing metadata columns
# - Session export to outputs/scRNA/sessioninfo/
# ==========================================================

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

mod_sc_qc_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    h3("Quality control"),
    br(),
    
    bs4Dash::bs4Card(
      title = "Help: how to use quality control",
      width = 12,
      solidHeader = TRUE,
      status = "info",
      collapsible = TRUE,
      collapsed = TRUE,
      icon = icon("circle-info"),
      tags$h4("Purpose"),
      tags$p("This module helps remove low-quality cells, empty droplets, stressed cells, and possible doublets before normalization, PCA, UMAP, clustering, and marker detection."),
      tags$h4("Main QC metrics"),
      tags$ul(
        tags$li(tags$b("nFeature_RNA:"), " number of detected genes per cell. Very low values often indicate poor-quality cells. Very high values may indicate doublets."),
        tags$li(tags$b("nCount_RNA:"), " total RNA counts or UMIs per cell. Very low values often indicate empty droplets. Very high values may indicate doublets or multiplets."),
        tags$li(tags$b("percent.mito:"), " mitochondrial read percentage. High values may indicate damaged, stressed, or dying cells."),
        tags$li(tags$b("percent.ribo:"), " ribosomal read percentage. Strong shifts may indicate technical effects or sample-specific biology.")
      ),
      tags$h4("Recommended workflow"),
      tags$ol(
        tags$li("Inspect violin plots and scatter plots before filtering."),
        tags$li("Adjust sliders or numeric thresholds based on your dataset distribution."),
        tags$li("Use the retained/removed cell summary to check how many cells will be filtered."),
        tags$li("Apply QC only when the selected thresholds look reasonable."),
        tags$li("Continue to normalization and PCA after QC filtering.")
      ),
      tags$h4("Important caution"),
      tags$p("Do not use one universal threshold for every dataset. Thresholds depend on tissue, species, sequencing depth, cell type diversity, and sample quality.")
    ),
    
    fluidRow(
      bs4Dash::bs4Card(
        title = "Filter low-quality cells",
        width = 12,
        solidHeader = TRUE,
        status = "primary",
        
        fluidRow(
          column(3, uiOutput(ns("slider_nb_features"))),
          column(3, uiOutput(ns("slider_nb_UMI"))),
          column(3, uiOutput(ns("slider_percent_Mito"))),
          column(3, uiOutput(ns("slider_percent_Ribo")))
        ),
        br(),
        fluidRow(
          column(3, actionButton(ns("updateF"), "Update from sliders", class = "btn-primary")),
          column(3, actionButton(ns("updateF_num"), "Update from numeric inputs", class = "btn-secondary")),
          column(6, uiOutput(ns("qc_filter_summary")))
        ),
        br(),
        fluidRow(
          column(3, h5("nFeature_RNA"), uiOutput(ns("text_nb_featuresMin")), uiOutput(ns("text_nb_featuresMax"))),
          column(3, h5("nCount_RNA"), uiOutput(ns("text_nb_UMIMin")), uiOutput(ns("text_nb_UMIMax"))),
          column(3, h5("percent.mito"), uiOutput(ns("text_percent_MitoMin")), uiOutput(ns("text_percent_MitoMax"))),
          column(3, h5("percent.ribo"), uiOutput(ns("text_percent_RiboMin")), uiOutput(ns("text_percent_RiboMax")))
        ),
        br(),
        fluidRow(
          column(4, uiOutput(ns("button_QC"))),
          column(4, downloadButton(ns("download_qc_session"), "Save QC session (.rds)", class = "btn-success")),
          column(4, uiOutput(ns("qc_threshold_warning")))
        )
      )
    ),
    
    bs4Dash::bs4Card(
      title = "Interpretation: QC filtering thresholds",
      width = 12,
      solidHeader = TRUE,
      status = "info",
      collapsible = TRUE,
      collapsed = TRUE,
      icon = icon("circle-info"),
      tags$h4("How to choose thresholds"),
      tags$ul(
        tags$li("Start with broad thresholds, then narrow them after inspecting plots."),
        tags$li("Remove cells in extreme low nFeature_RNA and nCount_RNA tails."),
        tags$li("Remove cells with unusually high percent.mito if they form a clear stressed-cell tail."),
        tags$li("Be careful with high nFeature_RNA or nCount_RNA cutoffs because some large or active cells naturally contain more RNA."),
        tags$li("Use doublet detection later if high-count cells remain suspicious.")
      ),
      tags$h4("What the colors mean"),
      tags$p("Cells marked as retained pass the current thresholds. Cells marked as removed fail at least one threshold and will be excluded when you apply QC."),
      tags$h4("Recommended action"),
      tags$p("After applying QC, re-check the plots. The remaining cells should show fewer extreme outliers without removing expected biological populations.")
    ),
    
    fluidRow(
      bs4Dash::bs4Card(
        title = "QC visualization",
        width = 12,
        solidHeader = TRUE,
        status = "primary",
        
        tabsetPanel(
          tabPanel(
            "Violin plots",
            uiOutput(ns("select_factor")),
            br(),
            fluidRow(
              column(
                6,
                shinycssloaders::withSpinner(plotly::plotlyOutput(ns("Vln_nFeatures")), type = 4),
                br(),
                fluidRow(
                  column(6, downloadButton(ns("dl_vln_feat_pdf"), "Features PDF")),
                  column(6, downloadButton(ns("dl_vln_feat_svg"), "Features SVG"))
                )
              ),
              column(
                6,
                shinycssloaders::withSpinner(plotly::plotlyOutput(ns("Vln_nCount")), type = 4),
                br(),
                fluidRow(
                  column(6, downloadButton(ns("dl_vln_umi_pdf"), "UMI PDF")),
                  column(6, downloadButton(ns("dl_vln_umi_svg"), "UMI SVG"))
                )
              )
            ),
            br(),
            fluidRow(
              column(
                6,
                shinycssloaders::withSpinner(plotly::plotlyOutput(ns("Vln_percent_Mito")), type = 4),
                br(),
                fluidRow(
                  column(6, downloadButton(ns("dl_vln_mito_pdf"), "Mito PDF")),
                  column(6, downloadButton(ns("dl_vln_mito_svg"), "Mito SVG"))
                )
              ),
              column(
                6,
                shinycssloaders::withSpinner(plotly::plotlyOutput(ns("Vln_percent_Ribo")), type = 4),
                br(),
                fluidRow(
                  column(6, downloadButton(ns("dl_vln_ribo_pdf"), "Ribo PDF")),
                  column(6, downloadButton(ns("dl_vln_ribo_svg"), "Ribo SVG"))
                )
              )
            ),
            br(),
            bs4Dash::bs4Card(
              title = "Interpretation: violin plots",
              width = 12,
              solidHeader = TRUE,
              status = "info",
              collapsible = TRUE,
              collapsed = TRUE,
              icon = icon("circle-info"),
              tags$p("Violin plots show the distribution of each QC metric across cells or metadata groups."),
              tags$ul(
                tags$li("Wide regions show where many cells are concentrated."),
                tags$li("Long low tails in nFeature_RNA or nCount_RNA may indicate poor-quality cells."),
                tags$li("Long high tails in nFeature_RNA or nCount_RNA may indicate doublets or large transcriptionally active cells."),
                tags$li("High percent.mito tails may indicate stressed or damaged cells."),
                tags$li("Group-specific shifts may indicate sample quality differences or real biological differences.")
              )
            )
          ),
          tabPanel(
            "Scatter",
            fluidRow(
              column(6, uiOutput(ns("select_paramX"))),
              column(6, uiOutput(ns("select_paramY")))
            ),
            br(),
            shinycssloaders::withSpinner(plotly::plotlyOutput(ns("scaterplot")), type = 4),
            br(),
            fluidRow(
              column(6, downloadButton(ns("dl_scatter_pdf"), "Scatter PDF")),
              column(6, downloadButton(ns("dl_scatter_svg"), "Scatter SVG"))
            ),
            br(),
            bs4Dash::bs4Card(
              title = "Interpretation: scatter plots",
              width = 12,
              solidHeader = TRUE,
              status = "info",
              collapsible = TRUE,
              collapsed = TRUE,
              icon = icon("circle-info"),
              tags$p("Scatter plots help identify relationships between QC metrics and reveal outlier groups."),
              tags$ul(
                tags$li("nCount_RNA and nFeature_RNA usually increase together."),
                tags$li("Cells with high nCount_RNA but unusually low nFeature_RNA may need inspection."),
                tags$li("Cells with high percent.mito and low nFeature_RNA are often low quality."),
                tags$li("Separated outlier clouds may represent technical artifacts, damaged cells, or a distinct biological population."),
                tags$li("Use marker genes and metadata before removing a cluster that may represent real biology.")
              )
            )
          )
        )
      )
    )
  )
}

mod_sc_qc_server <- function(id, seurat_r, sc_state = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    current_seu <- reactiveVal(NULL)
    filtered_seu <- reactiveVal(NULL)
    
    params <- reactiveValues(
      min1 = NULL, max1 = NULL,
      min2 = NULL, max2 = NULL,
      min3 = NULL, max3 = NULL,
      min4 = NULL, max4 = NULL
    )
    
    safe_theme <- function() {
      if (exists("theme_cotra", mode = "function")) {
        theme_cotra()
      } else {
        ggplot2::theme_bw()
      }
    }
    
    clean_numeric <- function(x, fallback = NA_real_) {
      out <- suppressWarnings(as.numeric(x))
      if (length(out) == 0 || is.na(out)) return(fallback)
      out
    }
    
    metric_range <- function(vals, round_to = 100, lower_zero = FALSE) {
      vals <- vals[is.finite(vals)]
      if (length(vals) == 0) return(c(0, 1))
      mn <- min(vals, na.rm = TRUE)
      mx <- max(vals, na.rm = TRUE)
      if (lower_zero) {
        mn <- 0
        mx <- ceiling(mx)
      } else {
        mn <- floor(mn / round_to) * round_to
        mx <- ceiling(mx / round_to) * round_to
      }
      if (mn == mx) mx <- mx + 1
      c(mn, mx)
    }
    
    valid_thresholds <- reactive({
      out <- character(0)
      pairs <- list(
        nFeature_RNA = c(params$min1, params$max1),
        nCount_RNA = c(params$min2, params$max2),
        percent.mito = c(params$min3, params$max3),
        percent.ribo = c(params$min4, params$max4)
      )
      
      for (nm in names(pairs)) {
        vals <- suppressWarnings(as.numeric(pairs[[nm]]))
        if (length(vals) == 2 && all(!is.na(vals)) && vals[1] > vals[2]) {
          out <- c(out, paste0(nm, " minimum is larger than maximum."))
        }
      }
      out
    })
    
    observe({
      seu <- seurat_r()
      if (is.null(seu)) return(NULL)
      current_seu(seu)
      
      if (!is.null(sc_state)) {
        sc_state$seurat <- seu
      }
      
      md <- seu@meta.data
      
      if ("nFeature_RNA" %in% colnames(md)) {
        rng <- metric_range(md$nFeature_RNA, round_to = 100)
        if (is.null(params$min1)) params$min1 <- rng[1]
        if (is.null(params$max1)) params$max1 <- rng[2]
      }
      if ("nCount_RNA" %in% colnames(md)) {
        rng <- metric_range(md$nCount_RNA, round_to = 100)
        if (is.null(params$min2)) params$min2 <- rng[1]
        if (is.null(params$max2)) params$max2 <- rng[2]
      }
      if ("percent.mito" %in% colnames(md)) {
        rng <- metric_range(md$percent.mito, lower_zero = TRUE)
        if (is.null(params$min3)) params$min3 <- rng[1]
        if (is.null(params$max3)) params$max3 <- rng[2]
      }
      if ("percent.ribo" %in% colnames(md)) {
        rng <- metric_range(md$percent.ribo, lower_zero = TRUE)
        if (is.null(params$min4)) params$min4 <- rng[1]
        if (is.null(params$max4)) params$max4 <- rng[2]
      }
    })
    
    output$slider_nb_features <- renderUI({
      seu <- current_seu()
      if (is.null(seu)) return(NULL)
      md <- seu@meta.data
      if (!"nFeature_RNA" %in% colnames(md)) return(NULL)
      rng <- metric_range(md$nFeature_RNA, round_to = 100)
      sliderInput(ns("qc_nb_feature"), "Number of genes per cell", min = rng[1], max = rng[2], step = 1, value = c(params$min1 %||% rng[1], params$max1 %||% rng[2]))
    })
    
    output$text_nb_featuresMin <- renderUI({
      seu <- current_seu()
      if (is.null(seu)) return(NULL)
      md <- seu@meta.data
      if (!"nFeature_RNA" %in% colnames(md)) return(NULL)
      rng <- metric_range(md$nFeature_RNA, round_to = 100)
      textInput(ns("qc_nb_featMin"), "Min", value = params$min1 %||% rng[1])
    })
    
    output$text_nb_featuresMax <- renderUI({
      seu <- current_seu()
      if (is.null(seu)) return(NULL)
      md <- seu@meta.data
      if (!"nFeature_RNA" %in% colnames(md)) return(NULL)
      rng <- metric_range(md$nFeature_RNA, round_to = 100)
      textInput(ns("qc_nb_featMax"), "Max", value = params$max1 %||% rng[2])
    })
    
    output$slider_nb_UMI <- renderUI({
      seu <- current_seu()
      if (is.null(seu)) return(NULL)
      md <- seu@meta.data
      if (!"nCount_RNA" %in% colnames(md)) return(NULL)
      rng <- metric_range(md$nCount_RNA, round_to = 100)
      sliderInput(ns("qc_nb_UMI"), "Number of UMI per cell", min = rng[1], max = rng[2], step = 1, value = c(params$min2 %||% rng[1], params$max2 %||% rng[2]))
    })
    
    output$text_nb_UMIMin <- renderUI({
      seu <- current_seu()
      if (is.null(seu)) return(NULL)
      md <- seu@meta.data
      if (!"nCount_RNA" %in% colnames(md)) return(NULL)
      rng <- metric_range(md$nCount_RNA, round_to = 100)
      textInput(ns("qc_nb_UMIMin"), "Min", value = params$min2 %||% rng[1])
    })
    
    output$text_nb_UMIMax <- renderUI({
      seu <- current_seu()
      if (is.null(seu)) return(NULL)
      md <- seu@meta.data
      if (!"nCount_RNA" %in% colnames(md)) return(NULL)
      rng <- metric_range(md$nCount_RNA, round_to = 100)
      textInput(ns("qc_nb_UMIMax"), "Max", value = params$max2 %||% rng[2])
    })
    
    output$slider_percent_Mito <- renderUI({
      seu <- current_seu()
      if (is.null(seu)) return(NULL)
      md <- seu@meta.data
      if (!"percent.mito" %in% colnames(md)) return(NULL)
      rng <- metric_range(md$percent.mito, lower_zero = TRUE)
      sliderInput(ns("qc_percent_Mito"), "Percent mito per cell", min = rng[1], max = rng[2], step = 0.1, value = c(params$min3 %||% rng[1], params$max3 %||% rng[2]))
    })
    
    output$text_percent_MitoMin <- renderUI({
      seu <- current_seu()
      if (is.null(seu)) return(NULL)
      md <- seu@meta.data
      if (!"percent.mito" %in% colnames(md)) return(NULL)
      rng <- metric_range(md$percent.mito, lower_zero = TRUE)
      textInput(ns("qc_percent_MitoMin"), "Min", value = params$min3 %||% rng[1])
    })
    
    output$text_percent_MitoMax <- renderUI({
      seu <- current_seu()
      if (is.null(seu)) return(NULL)
      md <- seu@meta.data
      if (!"percent.mito" %in% colnames(md)) return(NULL)
      rng <- metric_range(md$percent.mito, lower_zero = TRUE)
      textInput(ns("qc_percent_MitoMax"), "Max", value = params$max3 %||% rng[2])
    })
    
    output$slider_percent_Ribo <- renderUI({
      seu <- current_seu()
      if (is.null(seu)) return(NULL)
      md <- seu@meta.data
      if (!"percent.ribo" %in% colnames(md)) return(NULL)
      rng <- metric_range(md$percent.ribo, lower_zero = TRUE)
      sliderInput(ns("qc_percent_Ribo"), "Percent ribo per cell", min = rng[1], max = rng[2], step = 0.1, value = c(params$min4 %||% rng[1], params$max4 %||% rng[2]))
    })
    
    output$text_percent_RiboMin <- renderUI({
      seu <- current_seu()
      if (is.null(seu)) return(NULL)
      md <- seu@meta.data
      if (!"percent.ribo" %in% colnames(md)) return(NULL)
      rng <- metric_range(md$percent.ribo, lower_zero = TRUE)
      textInput(ns("qc_percent_RiboMin"), "Min", value = params$min4 %||% rng[1])
    })
    
    output$text_percent_RiboMax <- renderUI({
      seu <- current_seu()
      if (is.null(seu)) return(NULL)
      md <- seu@meta.data
      if (!"percent.ribo" %in% colnames(md)) return(NULL)
      rng <- metric_range(md$percent.ribo, lower_zero = TRUE)
      textInput(ns("qc_percent_RiboMax"), "Max", value = params$max4 %||% rng[2])
    })
    
    factorCol <- reactive({
      seu <- current_seu()
      if (is.null(seu)) return(NULL)
      md <- seu@meta.data
      cols <- colnames(md)[vapply(md, function(x) is.factor(x) || is.character(x), logical(1))]
      if (length(cols) == 0) {
        return("All_cells")
      }
      cols
    })
    
    output$select_factor <- renderUI({
      choicesFactor <- factorCol()
      if (is.null(choicesFactor)) return(NULL)
      selectInput(ns("color_factor"), "Metadata for grouping", choices = choicesFactor, selected = choicesFactor[1], width = "40%")
    })
    
    df_all <- reactive({
      seu <- current_seu()
      if (is.null(seu)) return(NULL)
      md <- seu@meta.data
      
      color_factor <- input$color_factor
      if (is.null(color_factor) || color_factor == "All_cells" || !color_factor %in% colnames(md)) {
        ident <- factor(rep("All cells", nrow(md)))
      } else {
        ident <- as.factor(md[[color_factor]])
      }
      
      temp <- data.frame(
        Cells_to_conserved = TRUE,
        Ident = ident,
        row.names = rownames(md)
      )
      
      if ("nFeature_RNA" %in% colnames(md)) {
        vals <- md$nFeature_RNA
        rng <- metric_range(vals, round_to = 100)
        min1 <- clean_numeric(params$min1, rng[1])
        max1 <- clean_numeric(params$max1, rng[2])
        temp$nFeature_RNA <- vals
        rem <- which(vals < min1 | vals > max1)
        if (length(rem) > 0) temp$Cells_to_conserved[rem] <- FALSE
      }
      
      if ("nCount_RNA" %in% colnames(md)) {
        vals <- md$nCount_RNA
        rng <- metric_range(vals, round_to = 100)
        min2 <- clean_numeric(params$min2, rng[1])
        max2 <- clean_numeric(params$max2, rng[2])
        temp$nCount_RNA <- vals
        rem <- which(vals < min2 | vals > max2)
        if (length(rem) > 0) temp$Cells_to_conserved[rem] <- FALSE
      }
      
      if ("percent.mito" %in% colnames(md)) {
        vals <- md$percent.mito
        rng <- metric_range(vals, lower_zero = TRUE)
        min3 <- clean_numeric(params$min3, rng[1])
        max3 <- clean_numeric(params$max3, rng[2])
        temp$percent.mito <- vals
        rem <- which(vals < min3 | vals > max3)
        if (length(rem) > 0) temp$Cells_to_conserved[rem] <- FALSE
      }
      
      if ("percent.ribo" %in% colnames(md)) {
        vals <- md$percent.ribo
        rng <- metric_range(vals, lower_zero = TRUE)
        min4 <- clean_numeric(params$min4, rng[1])
        max4 <- clean_numeric(params$max4, rng[2])
        temp$percent.ribo <- vals
        rem <- which(vals < min4 | vals > max4)
        if (length(rem) > 0) temp$Cells_to_conserved[rem] <- FALSE
      }
      
      temp
    })
    
    output$qc_filter_summary <- renderUI({
      df <- df_all()
      if (is.null(df)) return(NULL)
      total <- nrow(df)
      kept <- sum(df$Cells_to_conserved)
      removed <- total - kept
      pct <- round((removed / total) * 100, 2)
      tags$div(
        style = "padding-top: 8px;",
        tags$b("Current filter preview: "),
        paste0(kept, " retained, ", removed, " removed, ", pct, "% removed")
      )
    })
    
    output$qc_threshold_warning <- renderUI({
      warnings <- valid_thresholds()
      if (length(warnings) == 0) return(NULL)
      tags$div(
        style = "color:#b00020; padding-top:8px;",
        tags$b("Threshold warning: "),
        paste(warnings, collapse = " ")
      )
    })
    
    make_violin <- function(df, ycol, ylab_text) {
      if (is.null(df) || !ycol %in% colnames(df)) return(NULL)
      ggplot2::ggplot(df, ggplot2::aes(x = Ident, y = .data[[ycol]])) +
        ggplot2::geom_violin(ggplot2::aes(fill = Ident), show.legend = FALSE) +
        ggplot2::geom_jitter(height = 0, size = 0.15, ggplot2::aes(color = Cells_to_conserved)) +
        ggplot2::scale_color_manual(values = c("FALSE" = "grey", "TRUE" = "black")) +
        ggplot2::xlab("Group") +
        ggplot2::ylab(ylab_text) +
        safe_theme()
    }
    
    qc_plot1_base <- reactive(make_violin(df_all(), "nFeature_RNA", "nFeature_RNA"))
    qc_plot2_base <- reactive(make_violin(df_all(), "nCount_RNA", "nCount_RNA"))
    qc_plot3_base <- reactive(make_violin(df_all(), "percent.mito", "percent.mito"))
    qc_plot4_base <- reactive(make_violin(df_all(), "percent.ribo", "percent.ribo"))
    
    output$Vln_nFeatures <- plotly::renderPlotly({
      p <- qc_plot1_base()
      if (is.null(p)) return(NULL)
      plotly::ggplotly(p)
    })
    
    output$Vln_nCount <- plotly::renderPlotly({
      p <- qc_plot2_base()
      if (is.null(p)) return(NULL)
      plotly::ggplotly(p)
    })
    
    output$Vln_percent_Mito <- plotly::renderPlotly({
      p <- qc_plot3_base()
      if (is.null(p)) return(NULL)
      plotly::ggplotly(p)
    })
    
    output$Vln_percent_Ribo <- plotly::renderPlotly({
      p <- qc_plot4_base()
      if (is.null(p)) return(NULL)
      plotly::ggplotly(p)
    })
    
    output$select_paramX <- renderUI({
      seu <- current_seu()
      if (is.null(seu)) return(NULL)
      md <- seu@meta.data
      choices <- c()
      if ("nFeature_RNA" %in% colnames(md)) choices <- c(choices, "Nb_Gene")
      if ("nCount_RNA" %in% colnames(md)) choices <- c(choices, "Nb_UMI")
      if ("percent.mito" %in% colnames(md)) choices <- c(choices, "Percent_mito")
      if ("percent.ribo" %in% colnames(md)) choices <- c(choices, "Percent_ribo")
      if (length(choices) == 0) return(NULL)
      selectInput(ns("paramX"), "X axis", choices = choices, selected = choices[1])
    })
    
    output$select_paramY <- renderUI({
      seu <- current_seu()
      if (is.null(seu)) return(NULL)
      md <- seu@meta.data
      choices <- c()
      if ("nFeature_RNA" %in% colnames(md)) choices <- c(choices, "Nb_Gene")
      if ("nCount_RNA" %in% colnames(md)) choices <- c(choices, "Nb_UMI")
      if ("percent.mito" %in% colnames(md)) choices <- c(choices, "Percent_mito")
      if ("percent.ribo" %in% colnames(md)) choices <- c(choices, "Percent_ribo")
      if (length(choices) == 0) return(NULL)
      selected <- if (length(choices) >= 2) choices[2] else choices[1]
      selectInput(ns("paramY"), "Y axis", choices = choices, selected = selected)
    })
    
    paramX <- reactive({
      if (is.null(input$paramX)) return(NULL)
      switch(input$paramX,
             "Nb_Gene" = "nFeature_RNA",
             "Nb_UMI" = "nCount_RNA",
             "Percent_mito" = "percent.mito",
             "Percent_ribo" = "percent.ribo",
             NULL
      )
    })
    
    paramY <- reactive({
      if (is.null(input$paramY)) return(NULL)
      switch(input$paramY,
             "Nb_Gene" = "nFeature_RNA",
             "Nb_UMI" = "nCount_RNA",
             "Percent_mito" = "percent.mito",
             "Percent_ribo" = "percent.ribo",
             NULL
      )
    })
    
    sca_plot_base <- reactive({
      df <- df_all()
      if (is.null(df)) return(NULL)
      xcol <- paramX()
      ycol <- paramY()
      if (is.null(xcol) || is.null(ycol)) return(NULL)
      if (!xcol %in% colnames(df) || !ycol %in% colnames(df)) return(NULL)
      
      ggplot2::ggplot(df, ggplot2::aes(x = .data[[xcol]], y = .data[[ycol]])) +
        ggplot2::geom_point(ggplot2::aes(color = Cells_to_conserved), size = 0.3) +
        ggplot2::scale_color_manual(values = c("FALSE" = "#D9717D", "TRUE" = "#4DB6D0")) +
        ggplot2::xlab(input$paramX) +
        ggplot2::ylab(input$paramY) +
        safe_theme()
    })
    
    output$scaterplot <- plotly::renderPlotly({
      p <- sca_plot_base()
      if (is.null(p)) return(NULL)
      plotly::ggplotly(p)
    })
    
    output$button_QC <- renderUI({
      df <- df_all()
      if (is.null(df)) return(NULL)
      if (length(valid_thresholds()) > 0) return(NULL)
      n_rm <- sum(!df$Cells_to_conserved)
      if (n_rm <= 0) return(NULL)
      actionButton(ns("apply_qc"), label = paste0("Filter out ", n_rm, " low-quality cells"), class = "btn-danger")
    })
    
    observeEvent(input$apply_qc, {
      seu <- current_seu()
      df <- df_all()
      if (is.null(seu) || is.null(df)) return(NULL)
      
      if (length(valid_thresholds()) > 0) {
        showModal(modalDialog(title = "Invalid thresholds", paste(valid_thresholds(), collapse = " "), easyClose = TRUE))
        return(NULL)
      }
      
      keep <- rownames(df)[df$Cells_to_conserved]
      if (length(keep) < 50) {
        showModal(modalDialog(title = "Too few cells", "Less than 50 cells would remain. Relax your thresholds.", easyClose = TRUE))
        return(NULL)
      }
      
      seu_qc <- subset(seu, cells = keep)
      seu_qc@meta.data$CoTRA_QC_pass <- TRUE
      seu_qc@misc$CoTRA_qc <- list(
        qc_applied = TRUE,
        thresholds = list(
          nFeature_RNA = c(params$min1, params$max1),
          nCount_RNA = c(params$min2, params$max2),
          percent.mito = c(params$min3, params$max3),
          percent.ribo = c(params$min4, params$max4)
        ),
        cells_before = ncol(seu),
        cells_after = ncol(seu_qc),
        cells_removed = ncol(seu) - ncol(seu_qc),
        created_at = Sys.time()
      )
      
      filtered_seu(seu_qc)
      current_seu(seu_qc)
      
      if (!is.null(sc_state)) {
        sc_state$seurat <- seu_qc
        sc_state$parameters$qc <- seu_qc@misc$CoTRA_qc
      }
      
      showModal(modalDialog(title = "QC applied", paste0("Kept ", length(keep), " cells after filtering."), easyClose = TRUE))
    })
    
    observeEvent(input$updateF_num, {
      params$min1 <- clean_numeric(input$qc_nb_featMin, params$min1)
      params$max1 <- clean_numeric(input$qc_nb_featMax, params$max1)
      params$min2 <- clean_numeric(input$qc_nb_UMIMin, params$min2)
      params$max2 <- clean_numeric(input$qc_nb_UMIMax, params$max2)
      params$min3 <- clean_numeric(input$qc_percent_MitoMin, params$min3)
      params$max3 <- clean_numeric(input$qc_percent_MitoMax, params$max3)
      params$min4 <- clean_numeric(input$qc_percent_RiboMin, params$min4)
      params$max4 <- clean_numeric(input$qc_percent_RiboMax, params$max4)
    }, ignoreInit = TRUE)
    
    observeEvent(input$updateF, {
      if (!is.null(input$qc_nb_feature)) {
        params$min1 <- input$qc_nb_feature[1]
        params$max1 <- input$qc_nb_feature[2]
      }
      if (!is.null(input$qc_nb_UMI)) {
        params$min2 <- input$qc_nb_UMI[1]
        params$max2 <- input$qc_nb_UMI[2]
      }
      if (!is.null(input$qc_percent_Mito)) {
        params$min3 <- input$qc_percent_Mito[1]
        params$max3 <- input$qc_percent_Mito[2]
      }
      if (!is.null(input$qc_percent_Ribo)) {
        params$min4 <- input$qc_percent_Ribo[1]
        params$max4 <- input$qc_percent_Ribo[2]
      }
    }, ignoreInit = TRUE)
    
    observeEvent({ list(params$min1, params$max1) }, {
      if (!is.null(params$min1) && !is.null(params$max1)) {
        updateSliderInput(session, "qc_nb_feature", value = c(params$min1, params$max1))
        updateTextInput(session, "qc_nb_featMin", value = params$min1)
        updateTextInput(session, "qc_nb_featMax", value = params$max1)
      }
    }, ignoreInit = TRUE)
    
    observeEvent({ list(params$min2, params$max2) }, {
      if (!is.null(params$min2) && !is.null(params$max2)) {
        updateSliderInput(session, "qc_nb_UMI", value = c(params$min2, params$max2))
        updateTextInput(session, "qc_nb_UMIMin", value = params$min2)
        updateTextInput(session, "qc_nb_UMIMax", value = params$max2)
      }
    }, ignoreInit = TRUE)
    
    observeEvent({ list(params$min3, params$max3) }, {
      if (!is.null(params$min3) && !is.null(params$max3)) {
        updateSliderInput(session, "qc_percent_Mito", value = c(params$min3, params$max3))
        updateTextInput(session, "qc_percent_MitoMin", value = params$min3)
        updateTextInput(session, "qc_percent_MitoMax", value = params$max3)
      }
    }, ignoreInit = TRUE)
    
    observeEvent({ list(params$min4, params$max4) }, {
      if (!is.null(params$min4) && !is.null(params$max4)) {
        updateSliderInput(session, "qc_percent_Ribo", value = c(params$min4, params$max4))
        updateTextInput(session, "qc_percent_RiboMin", value = params$min4)
        updateTextInput(session, "qc_percent_RiboMax", value = params$max4)
      }
    }, ignoreInit = TRUE)
    
    save_plot <- function(plot_reactive, file, device = "pdf") {
      p <- plot_reactive()
      if (is.null(p)) return(NULL)
      ggplot2::ggsave(file, plot = p, device = device, width = 7, height = 5)
    }
    
    output$dl_vln_feat_pdf <- downloadHandler(filename = function() "QC_violin_features.pdf", content = function(file) save_plot(qc_plot1_base, file, "pdf"))
    output$dl_vln_feat_svg <- downloadHandler(filename = function() "QC_violin_features.svg", content = function(file) save_plot(qc_plot1_base, file, "svg"))
    output$dl_vln_umi_pdf <- downloadHandler(filename = function() "QC_violin_UMI.pdf", content = function(file) save_plot(qc_plot2_base, file, "pdf"))
    output$dl_vln_umi_svg <- downloadHandler(filename = function() "QC_violin_UMI.svg", content = function(file) save_plot(qc_plot2_base, file, "svg"))
    output$dl_vln_mito_pdf <- downloadHandler(filename = function() "QC_violin_mito.pdf", content = function(file) save_plot(qc_plot3_base, file, "pdf"))
    output$dl_vln_mito_svg <- downloadHandler(filename = function() "QC_violin_mito.svg", content = function(file) save_plot(qc_plot3_base, file, "svg"))
    output$dl_vln_ribo_pdf <- downloadHandler(filename = function() "QC_violin_ribo.pdf", content = function(file) save_plot(qc_plot4_base, file, "pdf"))
    output$dl_vln_ribo_svg <- downloadHandler(filename = function() "QC_violin_ribo.svg", content = function(file) save_plot(qc_plot4_base, file, "svg"))
    output$dl_scatter_pdf <- downloadHandler(filename = function() "QC_scatter.pdf", content = function(file) save_plot(sca_plot_base, file, "pdf"))
    output$dl_scatter_svg <- downloadHandler(filename = function() "QC_scatter.svg", content = function(file) save_plot(sca_plot_base, file, "svg"))
    
    output$download_qc_session <- downloadHandler(
      filename = function() {
        paste0("CoTRA_scRNA_QC_session_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds")
      },
      content = function(file) {
        out_dir <- file.path("outputs", "scRNA", "sessioninfo")
        if (!dir.exists(out_dir)) {
          dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
        }
        
        out_file <- file.path(
          out_dir,
          paste0("CoTRA_scRNA_QC_session_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds")
        )
        
        seu_out <- if (!is.null(filtered_seu())) filtered_seu() else current_seu()
        req(seu_out)
        
        qc_info <- list(
          step = "Quality_Control",
          saved_at = Sys.time(),
          seurat = seu_out,
          qc_applied = !is.null(filtered_seu()),
          qc_thresholds = list(
            nFeature_RNA = c(params$min1, params$max1),
            nCount_RNA = c(params$min2, params$max2),
            percent.mito = c(params$min3, params$max3),
            percent.ribo = c(params$min4, params$max4)
          ),
          cell_summary = list(
            total_cells_current = ncol(seu_out),
            total_genes_current = nrow(seu_out)
          ),
          filter_preview = df_all(),
          cotra_module = "modules/sc/qc.R",
          output_path = out_file
        )
        
        saveRDS(qc_info, out_file)
        file.copy(out_file, file, overwrite = TRUE)
      }
    )
    
    list(
      seurat = reactive({
        if (!is.null(filtered_seu())) filtered_seu() else seurat_r()
      }),
      qc_applied = reactive(!is.null(filtered_seu())),
      qc_thresholds = reactive({
        list(
          nFeature_RNA = c(params$min1, params$max1),
          nCount_RNA = c(params$min2, params$max2),
          percent.mito = c(params$min3, params$max3),
          percent.ribo = c(params$min4, params$max4)
        )
      }),
      session_path = reactive(file.path("outputs", "scRNA", "sessioninfo"))
    )
  })
}

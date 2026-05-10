# ==========================================================
# modules/sc/pca.R
# CoTRA scRNA-seq PCA Module
#
# Features:
# - Help and interpretation sections
# - Seurat v4/v5 safe scaling and PCA
# - Uses selected genes from gene_selection.R when available
# - Optional fallback to VariableFeatures
# - Elbow plot
# - PCA scatter plot
# - PC loading table
# - Top loading genes heatmap
# - Recommended number of PCs
# - Downloadable plots and tables
# - Downstream-safe return object for UMAP/t-SNE/clustering
# ==========================================================

mod_sc_pca_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    h3("Principal component analysis"),
    br(),
    
    bs4Dash::bs4Card(
      title = "Help: how to use PCA",
      width = 12,
      solidHeader = TRUE,
      status = "info",
      collapsible = TRUE,
      collapsed = TRUE,
      icon = icon("circle-info"),
      tags$h4("Purpose"),
      tags$p("PCA reduces thousands of genes into a smaller set of principal components. These components capture the strongest expression patterns across cells and are used for UMAP, t-SNE, clustering, and downstream analysis."),
      tags$h4("Recommended workflow"),
      tags$ol(
        tags$li("Run QC first."),
        tags$li("Run gene selection to select HVGs or a custom gene set."),
        tags$li("Run PCA using selected genes."),
        tags$li("Inspect the elbow plot and loading table."),
        tags$li("Choose a reasonable number of PCs for UMAP, t-SNE, and clustering.")
      ),
      tags$h4("Input genes"),
      tags$ul(
        tags$li(tags$b("Selected genes:"), " genes saved by the gene_selection.R module."),
        tags$li(tags$b("VariableFeatures fallback:"), " uses Seurat VariableFeatures if selected genes are not available."),
        tags$li(tags$b("All genes fallback:"), " can be used if no selected or variable genes exist, but this is slower and less recommended.")
      ),
      tags$h4("Important caution"),
      tags$p("PCA reflects the gene set used. If you use only custom genes, PCA will focus only on that biology and may miss other cell populations or technical effects.")
    ),
    
    fluidRow(
      bs4Dash::bs4Card(
        title = "PCA settings",
        width = 12,
        solidHeader = TRUE,
        status = "primary",
        fluidRow(
          column(
            3,
            selectInput(
              ns("assay_use"),
              "Assay",
              choices = c("RNA"),
              selected = "RNA"
            )
          ),
          column(
            3,
            selectInput(
              ns("gene_source"),
              "Genes for PCA",
              choices = c(
                "Use selected genes from gene selection" = "selected",
                "Use Seurat VariableFeatures" = "variable",
                "Use all genes" = "all"
              ),
              selected = "selected"
            )
          ),
          column(
            3,
            numericInput(
              ns("npcs"),
              "Number of PCs to compute",
              value = 50,
              min = 5,
              max = 100,
              step = 5
            )
          ),
          column(
            3,
            numericInput(
              ns("top_loading_genes"),
              "Top loading genes per PC",
              value = 15,
              min = 5,
              max = 50,
              step = 5
            )
          )
        ),
        fluidRow(
          column(3, checkboxInput(ns("regress_mito"), "Regress percent.mito", FALSE)),
          column(3, checkboxInput(ns("regress_ribo"), "Regress percent.ribo", FALSE)),
          column(3, checkboxInput(ns("regress_counts"), "Regress nCount_RNA", FALSE)),
          column(3, br(), actionButton(ns("run_pca"), "Run PCA", class = "btn-primary"))
        ),
        hr(),
        fluidRow(
          column(4, downloadButton(ns("download_pca_session"), "Save PCA session (.rds)", class = "btn-success")),
          column(8, uiOutput(ns("pca_session_save_hint")))
        ),
        br(),
        uiOutput(ns("pca_status_box"))
      )
    ),
    
    bs4Dash::bs4Card(
      title = "Interpretation: PCA settings",
      width = 12,
      solidHeader = TRUE,
      status = "info",
      collapsible = TRUE,
      collapsed = TRUE,
      icon = icon("circle-info"),
      tags$h4("How to choose genes"),
      tags$ul(
        tags$li("Use selected HVGs for standard exploratory scRNA-seq analysis."),
        tags$li("Use VariableFeatures if the gene selection module was skipped."),
        tags$li("Use all genes only when needed because it can increase noise and runtime."),
        tags$li("Use custom selected genes only when you want PCA to focus on a specific process or marker panel.")
      ),
      tags$h4("Regression options"),
      tags$p("Regression can reduce technical effects before PCA. Use it carefully. Regressing true biological signals can weaken meaningful separation."),
      tags$ul(
        tags$li(tags$b("percent.mito:"), " useful if mitochondrial content is a technical stress signal."),
        tags$li(tags$b("percent.ribo:"), " useful if ribosomal content dominates technical variation."),
        tags$li(tags$b("nCount_RNA:"), " useful if sequencing depth strongly drives PC1 or PC2.")
      )
    ),
    
    fluidRow(
      bs4Dash::bs4Card(
        title = "PCA results",
        width = 12,
        solidHeader = TRUE,
        status = "primary",
        tabsetPanel(
          tabPanel(
            "Elbow plot",
            br(),
            uiOutput(ns("recommended_pc_box")),
            shinycssloaders::withSpinner(plotly::plotlyOutput(ns("elbow_plot")), type = 4),
            br(),
            fluidRow(
              column(4, downloadButton(ns("dl_elbow_pdf"), "Elbow PDF")),
              column(4, downloadButton(ns("dl_elbow_svg"), "Elbow SVG")),
              column(4, downloadButton(ns("dl_pc_variance"), "PC variance CSV"))
            ),
            br(),
            bs4Dash::bs4Card(
              title = "Interpretation: elbow plot",
              width = 12,
              solidHeader = TRUE,
              status = "info",
              collapsible = TRUE,
              collapsed = TRUE,
              icon = icon("circle-info"),
              tags$p("The elbow plot shows how much variation each PC explains. The useful number of PCs is often near the point where the curve starts to flatten."),
              tags$ul(
                tags$li("A sharp drop followed by a flat curve suggests fewer PCs may be enough."),
                tags$li("A gradual curve suggests the dataset may need more PCs."),
                tags$li("Use the recommended PC number as a starting point, not a fixed rule."),
                tags$li("Check UMAP, clustering, and marker results after choosing PCs.")
              )
            )
          ),
          tabPanel(
            "PCA scatter",
            br(),
            fluidRow(
              column(3, numericInput(ns("pc_x"), "PC on X axis", value = 1, min = 1, max = 100, step = 1)),
              column(3, numericInput(ns("pc_y"), "PC on Y axis", value = 2, min = 1, max = 100, step = 1)),
              column(6, uiOutput(ns("color_by_ui")))
            ),
            shinycssloaders::withSpinner(plotly::plotlyOutput(ns("pca_scatter")), type = 4),
            br(),
            fluidRow(
              column(6, downloadButton(ns("dl_pca_scatter_pdf"), "Scatter PDF")),
              column(6, downloadButton(ns("dl_pca_scatter_svg"), "Scatter SVG"))
            ),
            br(),
            bs4Dash::bs4Card(
              title = "Interpretation: PCA scatter",
              width = 12,
              solidHeader = TRUE,
              status = "info",
              collapsible = TRUE,
              collapsed = TRUE,
              icon = icon("circle-info"),
              tags$p("The PCA scatter plot shows cells positioned by selected PCs. Separation can reflect cell types, cell states, batch effects, sample differences, or technical variation."),
              tags$ul(
                tags$li("Color by sample, batch, condition, or QC metrics to identify technical drivers."),
                tags$li("Strong separation by condition may reflect biology or batch, depending on study design."),
                tags$li("Strong separation by percent.mito or nCount_RNA suggests technical variation may influence PCA."),
                tags$li("Do not assign cell types from PCA alone. Use markers and clustering.")
              )
            )
          ),
          tabPanel(
            "PC loadings",
            br(),
            fluidRow(
              column(4, numericInput(ns("loading_pc"), "PC for loading table", value = 1, min = 1, max = 100, step = 1)),
              column(4, downloadButton(ns("dl_loadings"), "Download loadings CSV"))
            ),
            br(),
            DT::DTOutput(ns("loading_table")),
            br(),
            bs4Dash::bs4Card(
              title = "Interpretation: PC loadings",
              width = 12,
              solidHeader = TRUE,
              status = "info",
              collapsible = TRUE,
              collapsed = TRUE,
              icon = icon("circle-info"),
              tags$p("PC loadings show which genes contribute most strongly to each principal component."),
              tags$ul(
                tags$li("Positive and negative loading genes define opposite directions of the same PC."),
                tags$li("Known marker genes can help interpret which cell type or state drives a PC."),
                tags$li("Mitochondrial, ribosomal, or stress genes among top loadings may indicate technical effects."),
                tags$li("Use PC loadings with UMAP, clustering, and marker genes for stronger interpretation.")
              )
            )
          ),
          tabPanel(
            "Loading heatmap",
            br(),
            fluidRow(
              column(4, numericInput(ns("heatmap_npcs"), "Number of PCs in heatmap", value = 6, min = 2, max = 20, step = 1)),
              column(4, downloadButton(ns("dl_heatmap_pdf"), "Heatmap PDF")),
              column(4, downloadButton(ns("dl_heatmap_svg"), "Heatmap SVG"))
            ),
            br(),
            shinycssloaders::withSpinner(plotOutput(ns("loading_heatmap"), height = "700px"), type = 4),
            br(),
            bs4Dash::bs4Card(
              title = "Interpretation: loading heatmap",
              width = 12,
              solidHeader = TRUE,
              status = "info",
              collapsible = TRUE,
              collapsed = TRUE,
              icon = icon("circle-info"),
              tags$p("The loading heatmap displays top genes that drive selected PCs. It helps identify whether PCs are driven by biological markers or technical gene groups."),
              tags$ul(
                tags$li("Distinct marker genes suggest biological structure."),
                tags$li("Repeated mitochondrial or ribosomal genes suggest technical influence."),
                tags$li("If many PCs are dominated by technical genes, revisit QC, regression, or gene selection.")
              )
            )
          )
        )
      )
    )
  )
}

mod_sc_pca_server <- function(id, seurat_r, sc_state = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    `%||%` <- function(x, y) {
      if (is.null(x)) y else x
    }
    
    pca_seu <- reactiveVal(NULL)
    pca_ready <- reactiveVal(FALSE)
    pca_features <- reactiveVal(character(0))
    pc_variance_df <- reactiveVal(NULL)
    recommended_pcs <- reactiveVal(NULL)
    elbow_plot_gg <- reactiveVal(NULL)
    scatter_plot_gg <- reactiveVal(NULL)
    heatmap_plot_gg <- reactiveVal(NULL)
    last_session_path <- reactiveVal(NULL)
    
    get_seu <- reactive({
      if (!is.null(sc_state) && !is.null(sc_state$seurat)) return(sc_state$seurat)
      seurat_r()
    })
    
    safe_theme <- function() {
      if (exists("theme_cotra", mode = "function")) {
        theme_cotra()
      } else {
        ggplot2::theme_bw()
      }
    }
    
    assay_choices <- reactive({
      seu <- get_seu()
      if (is.null(seu)) return(c("RNA"))
      assays <- tryCatch(Seurat::Assays(seu), error = function(e) "RNA")
      if (length(assays) == 0) assays <- "RNA"
      assays
    })
    
    observe({
      choices <- assay_choices()
      selected <- if ("SCT" %in% choices) "SCT" else if ("RNA" %in% choices) "RNA" else choices[1]
      updateSelectInput(session, "assay_use", choices = choices, selected = selected)
    })
    
    assay_exists <- function(seu, assay) {
      isTRUE(assay %in% Seurat::Assays(seu))
    }
    
    get_assay_matrix_safe <- function(seu, assay = "RNA", preferred = c("data", "counts", "scale.data")) {
      preferred <- unique(preferred)
      
      for (layer_name in preferred) {
        mat <- tryCatch(
          Seurat::GetAssayData(seu, assay = assay, layer = layer_name),
          error = function(e) NULL
        )
        if (!is.null(mat) && nrow(mat) > 0 && ncol(mat) > 0) {
          return(list(matrix = mat, source = paste0("layer:", layer_name)))
        }
      }
      
      for (slot_name in preferred) {
        mat <- tryCatch(
          Seurat::GetAssayData(seu, assay = assay, slot = slot_name),
          error = function(e) NULL
        )
        if (!is.null(mat) && nrow(mat) > 0 && ncol(mat) > 0) {
          return(list(matrix = mat, source = paste0("slot:", slot_name)))
        }
      }
      
      list(matrix = NULL, source = NA_character_)
    }
    
    assay_has_data <- function(seu, assay) {
      res <- get_assay_matrix_safe(seu, assay = assay, preferred = c("data"))
      !is.null(res$matrix)
    }
    
    get_selected_features <- function(seu, assay, gene_source) {
      available <- rownames(seu[[assay]])
      
      if (gene_source == "selected") {
        selected <- character(0)
        selected <- tryCatch(seu@misc$CoTRA_gene_selection$selected_genes, error = function(e) character(0))
        if (is.null(selected)) selected <- character(0)
        selected <- selected[selected %in% available]
        if (length(selected) > 0) return(unique(selected))
        
        vf <- tryCatch(Seurat::VariableFeatures(seu, assay = assay), error = function(e) character(0))
        vf <- vf[vf %in% available]
        if (length(vf) > 0) return(unique(vf))
      }
      
      if (gene_source == "variable") {
        vf <- tryCatch(Seurat::VariableFeatures(seu, assay = assay), error = function(e) character(0))
        vf <- vf[vf %in% available]
        if (length(vf) > 0) return(unique(vf))
      }
      
      unique(available)
    }
    
    get_pc_embeddings <- function(seu) {
      emb <- tryCatch(Seurat::Embeddings(seu, reduction = "pca"), error = function(e) NULL)
      emb
    }
    
    get_pc_loadings <- function(seu) {
      load <- tryCatch(Seurat::Loadings(seu, reduction = "pca"), error = function(e) NULL)
      load
    }
    
    calculate_pc_variance <- function(seu) {
      stdev <- tryCatch(seu[["pca"]]@stdev, error = function(e) NULL)
      if (is.null(stdev) || length(stdev) == 0) return(NULL)
      variance <- stdev^2
      pct <- variance / sum(variance) * 100
      cum <- cumsum(pct)
      data.frame(
        PC = seq_along(pct),
        percent_variance = pct,
        cumulative_variance = cum,
        stringsAsFactors = FALSE
      )
    }
    
    recommend_pcs <- function(var_df) {
      if (is.null(var_df) || nrow(var_df) < 3) return(10)
      pct <- var_df$percent_variance
      diff1 <- abs(diff(pct))
      elbow <- which(diff1 < 0.1)[1]
      if (is.na(elbow)) elbow <- which(var_df$cumulative_variance >= 80)[1]
      if (is.na(elbow)) elbow <- min(30, nrow(var_df))
      elbow <- max(5, min(elbow, nrow(var_df)))
      elbow
    }
    
    make_elbow_plot <- function(var_df, recommended) {
      if (is.null(var_df)) return(NULL)
      ggplot2::ggplot(var_df, ggplot2::aes(x = PC, y = percent_variance)) +
        ggplot2::geom_point(size = 2) +
        ggplot2::geom_line() +
        ggplot2::geom_vline(xintercept = recommended, linetype = "dashed") +
        ggplot2::labs(
          title = paste0("Elbow plot, recommended PCs: ", recommended),
          x = "Principal component",
          y = "Percent variance explained"
        ) +
        safe_theme()
    }
    
    make_scatter_plot <- reactive({
      seu <- pca_seu()
      if (is.null(seu)) return(NULL)
      emb <- get_pc_embeddings(seu)
      if (is.null(emb)) return(NULL)
      
      pcx <- as.integer(input$pc_x %||% 1)
      pcy <- as.integer(input$pc_y %||% 2)
      pcx <- max(1, min(pcx, ncol(emb)))
      pcy <- max(1, min(pcy, ncol(emb)))
      
      df <- as.data.frame(emb[, c(pcx, pcy), drop = FALSE])
      colnames(df) <- c("PCX", "PCY")
      df$cell <- rownames(df)
      
      md <- seu@meta.data
      color_by <- input$color_by %||% "none"
      if (!is.null(color_by) && color_by != "none" && color_by %in% colnames(md)) {
        df$color_value <- md[rownames(df), color_by]
        color_label <- color_by
      } else {
        df$color_value <- "Cells"
        color_label <- "Cells"
      }
      
      ggplot2::ggplot(df, ggplot2::aes(x = PCX, y = PCY, color = color_value, text = cell)) +
        ggplot2::geom_point(size = 0.6, alpha = 0.8) +
        ggplot2::labs(
          title = paste0("PCA scatter: PC", pcx, " vs PC", pcy),
          x = paste0("PC", pcx),
          y = paste0("PC", pcy),
          color = color_label
        ) +
        safe_theme()
    })
    
    output$color_by_ui <- renderUI({
      seu <- pca_seu()
      if (is.null(seu)) return(NULL)
      md <- seu@meta.data
      choices <- c("None" = "none", colnames(md))
      selected <- if ("orig.ident" %in% colnames(md)) "orig.ident" else "none"
      selectInput(ns("color_by"), "Color cells by", choices = choices, selected = selected)
    })
    
    output$pca_session_save_hint <- renderUI({
      if (!isTRUE(pca_ready()) || is.null(pca_seu())) {
        return(tags$span(style = "color:#777;", "Run PCA before saving this step."))
      }
      if (is.null(last_session_path())) {
        return(tags$span(style = "color:#777;", "Save the current PCA state as a CoTRA project RDS."))
      }
      tags$span(style = "color:#2E7D32;", paste("Last saved:", last_session_path()))
    })
    
    build_pca_project <- function(seu_out, out_file = NULL) {
      list(
        seurat = seu_out,
        project_meta = list(
          project_type = "CoTRA_scRNA_project",
          step = "PCA",
          saved_at = Sys.time(),
          cotra_module = "modules/sc/pca.R",
          output_path = out_file
        ),
        cells_meta = seu_out@meta.data,
        step = "PCA",
        saved_at = Sys.time(),
        genes_list = list(
          selected_genes = tryCatch(seu_out@misc$CoTRA_gene_selection$selected_genes, error = function(e) character(0)),
          hvgs = tryCatch(seu_out@misc$CoTRA_gene_selection$hvg_genes, error = function(e) character(0)),
          pca_features = pca_features()
        ),
        parameters = list(
          qc = tryCatch(seu_out@misc$CoTRA_qc, error = function(e) NULL),
          gene_selection = tryCatch(seu_out@misc$CoTRA_gene_selection, error = function(e) NULL),
          pca = tryCatch(seu_out@misc$CoTRA_pca, error = function(e) NULL)
        ),
        tables = list(
          pc_variance = pc_variance_df(),
          current_pc_loadings = tryCatch(loading_table_df(), error = function(e) NULL)
        ),
        cotra_version = if (exists("cotra_version")) cotra_version else "unknown"
      )
    }
    
    output$pca_status_box <- renderUI({
      seu <- pca_seu()
      feats <- pca_features()
      
      if (is.null(seu) || !isTRUE(pca_ready())) {
        return(tags$div(style = "color:#777;", "PCA has not been run yet."))
      }
      
      var_df <- pc_variance_df()
      rec <- recommended_pcs()
      tags$div(
        tags$b("PCA status: "), "completed", tags$br(),
        tags$b("Genes used: "), length(feats), tags$br(),
        tags$b("PCs computed: "), ifelse(is.null(var_df), NA, nrow(var_df)), tags$br(),
        tags$b("Recommended PCs: "), rec
      )
    })
    
    output$recommended_pc_box <- renderUI({
      rec <- recommended_pcs()
      var_df <- pc_variance_df()
      if (is.null(rec) || is.null(var_df)) {
        return(tags$div(style = "color:#777;", "Run PCA to calculate recommended PCs."))
      }
      cum <- round(var_df$cumulative_variance[rec], 2)
      tags$div(
        style = "padding:8px; background:#f7f7f7; border-radius:6px; margin-bottom:10px;",
        tags$b("Recommended number of PCs: "), rec,
        tags$br(),
        tags$span(paste0("Cumulative variance at PC", rec, ": ", cum, "%"))
      )
    })
    
    observeEvent(input$run_pca, {
      seu <- get_seu()
      req(seu)
      
      withProgress(message = "Running PCA", value = 0.1, {
        assay <- input$assay_use %||% "RNA"
        
        if (!assay_exists(seu, assay)) {
          showModal(modalDialog(
            title = "Assay not found",
            paste("Assay", assay, "was not found in the Seurat object."),
            easyClose = TRUE
          ))
          return(NULL)
        }
        
        Seurat::DefaultAssay(seu) <- assay
        
        if (!assay_has_data(seu, assay) && assay != "SCT") {
          incProgress(0.2, detail = "Normalizing data")
          seu <- Seurat::NormalizeData(
            seu,
            assay = assay,
            normalization.method = "LogNormalize",
            scale.factor = 10000,
            verbose = FALSE
          )
        }
        
        incProgress(0.35, detail = "Selecting PCA features")
        
        features <- get_selected_features(seu, assay, input$gene_source %||% "selected")
        features <- unique(features)
        features <- features[features %in% rownames(seu[[assay]])]
        
        if (length(features) < 50) {
          showModal(modalDialog(
            title = "Too few genes for PCA",
            paste("Only", length(features), "genes are available for PCA. Run gene selection again or use VariableFeatures/all genes."),
            easyClose = TRUE
          ))
          return(NULL)
        }
        
        npcs <- as.integer(input$npcs %||% 50)
        max_pcs <- min(npcs, ncol(seu) - 1, length(features) - 1)
        if (max_pcs < 2) {
          showModal(modalDialog(title = "PCA cannot run", "Too few cells or genes for PCA.", easyClose = TRUE))
          return(NULL)
        }
        
        vars_to_regress <- character(0)
        md <- seu@meta.data
        
        requested_regress <- c(
          if (isTRUE(input$regress_mito)) "percent.mito" else NULL,
          if (isTRUE(input$regress_ribo)) "percent.ribo" else NULL,
          if (isTRUE(input$regress_counts)) "nCount_RNA" else NULL
        )
        
        vars_to_regress <- requested_regress[requested_regress %in% colnames(md)]
        vars_to_regress <- vars_to_regress[vapply(vars_to_regress, function(v) {
          vals <- md[[v]]
          is.numeric(vals) && any(is.finite(vals)) && stats::sd(vals, na.rm = TRUE) > 0
        }, logical(1))]
        
        missing_regress <- setdiff(requested_regress, vars_to_regress)
        if (length(missing_regress) > 0) {
          showNotification(
            paste("Skipped regression variables not present or not variable:", paste(missing_regress, collapse = ", ")),
            type = "warning",
            duration = 6
          )
        }
        
        incProgress(0.55, detail = "Scaling data")
        
        if (assay == "SCT") {
          seu <- tryCatch(
            Seurat::ScaleData(seu, assay = assay, features = features, verbose = FALSE),
            error = function(e) {
              showModal(modalDialog(title = "ScaleData failed", conditionMessage(e), easyClose = TRUE))
              NULL
            }
          )
        } else {
          if (length(vars_to_regress) > 0) {
            seu <- tryCatch(
              Seurat::ScaleData(
                seu,
                assay = assay,
                features = features,
                vars.to.regress = vars_to_regress,
                verbose = FALSE
              ),
              error = function(e) {
                showModal(modalDialog(title = "ScaleData failed", conditionMessage(e), easyClose = TRUE))
                NULL
              }
            )
          } else {
            seu <- tryCatch(
              Seurat::ScaleData(
                seu,
                assay = assay,
                features = features,
                verbose = FALSE
              ),
              error = function(e) {
                showModal(modalDialog(title = "ScaleData failed", conditionMessage(e), easyClose = TRUE))
                NULL
              }
            )
          }
        }
        
        if (is.null(seu)) return(NULL)
        
        incProgress(0.75, detail = "Computing PCA")
        
        seu <- tryCatch(
          Seurat::RunPCA(
            seu,
            assay = assay,
            features = features,
            npcs = max_pcs,
            verbose = FALSE
          ),
          error = function(e) {
            showModal(modalDialog(title = "RunPCA failed", conditionMessage(e), easyClose = TRUE))
            NULL
          }
        )
        
        if (is.null(seu)) return(NULL)
        
        var_df <- calculate_pc_variance(seu)
        rec <- recommend_pcs(var_df)
        elbow <- make_elbow_plot(var_df, rec)
        
        seu@misc$CoTRA_pca <- list(
          assay = assay,
          features = features,
          n_features = length(features),
          npcs = max_pcs,
          recommended_pcs = rec,
          pc_variance = var_df,
          vars_to_regress = vars_to_regress,
          created_at = Sys.time()
        )
        
        pca_seu(seu)
        pca_ready(TRUE)
        pca_features(features)
        pc_variance_df(var_df)
        recommended_pcs(rec)
        elbow_plot_gg(elbow)
        
        if (!is.null(sc_state)) {
          sc_state$seurat <- seu
          sc_state$parameters$pca <- seu@misc$CoTRA_pca
        }
        
        incProgress(1)
      })
    })
    
    output$elbow_plot <- plotly::renderPlotly({
      p <- elbow_plot_gg()
      if (is.null(p)) return(NULL)
      plotly::ggplotly(p)
    })
    
    output$pca_scatter <- plotly::renderPlotly({
      p <- make_scatter_plot()
      scatter_plot_gg(p)
      if (is.null(p)) return(NULL)
      plotly::ggplotly(p, tooltip = c("text", "x", "y", "color"))
    })
    
    output$loading_table <- DT::renderDT({
      seu <- pca_seu()
      req(seu)
      load <- get_pc_loadings(seu)
      req(load)
      
      pc <- as.integer(input$loading_pc %||% 1)
      pc <- max(1, min(pc, ncol(load)))
      vals <- load[, pc]
      df <- data.frame(
        gene = names(vals),
        loading = as.numeric(vals),
        abs_loading = abs(as.numeric(vals)),
        direction = ifelse(vals >= 0, "positive", "negative"),
        PC = paste0("PC", pc),
        stringsAsFactors = FALSE
      )
      df <- df[order(-df$abs_loading), , drop = FALSE]
      DT::datatable(df, options = list(pageLength = 25, scrollX = TRUE))
    })
    
    loading_table_df <- reactive({
      seu <- pca_seu()
      req(seu)
      load <- get_pc_loadings(seu)
      req(load)
      pc <- as.integer(input$loading_pc %||% 1)
      pc <- max(1, min(pc, ncol(load)))
      vals <- load[, pc]
      df <- data.frame(
        gene = names(vals),
        loading = as.numeric(vals),
        abs_loading = abs(as.numeric(vals)),
        direction = ifelse(vals >= 0, "positive", "negative"),
        PC = paste0("PC", pc),
        stringsAsFactors = FALSE
      )
      df[order(-df$abs_loading), , drop = FALSE]
    })
    
    make_loading_heatmap <- reactive({
      seu <- pca_seu()
      if (is.null(seu)) return(NULL)
      load <- get_pc_loadings(seu)
      if (is.null(load)) return(NULL)
      
      n_pc <- as.integer(input$heatmap_npcs %||% 6)
      n_pc <- max(1, min(n_pc, ncol(load)))
      top_n <- as.integer(input$top_loading_genes %||% 15)
      top_n <- max(5, min(top_n, nrow(load)))
      
      genes <- unique(unlist(lapply(seq_len(n_pc), function(i) {
        vals <- abs(load[, i])
        names(sort(vals, decreasing = TRUE))[seq_len(min(top_n, length(vals)))]
      })))
      
      genes <- genes[genes %in% rownames(load)]
      if (length(genes) == 0) return(NULL)
      
      df <- as.data.frame(load[genes, seq_len(n_pc), drop = FALSE])
      df$gene <- rownames(df)
      long <- tidyr::pivot_longer(df, cols = -gene, names_to = "PC", values_to = "loading")
      long$gene <- factor(long$gene, levels = rev(unique(long$gene)))
      
      ggplot2::ggplot(long, ggplot2::aes(x = PC, y = gene, fill = loading)) +
        ggplot2::geom_tile() +
        ggplot2::labs(
          title = "Top PCA loading genes",
          x = "Principal component",
          y = "Gene"
        ) +
        safe_theme()
    })
    
    output$loading_heatmap <- renderPlot({
      p <- make_loading_heatmap()
      heatmap_plot_gg(p)
      if (is.null(p)) return(NULL)
      print(p)
    })
    
    save_plot <- function(plot_reactive, file, device = "pdf", width = 7, height = 5) {
      p <- plot_reactive()
      req(p)
      ggplot2::ggsave(file, plot = p, device = device, width = width, height = height)
    }
    
    output$dl_elbow_pdf <- downloadHandler(
      filename = function() paste0("CoTRA_PCA_elbow_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"),
      content = function(file) save_plot(elbow_plot_gg, file, "pdf")
    )
    
    output$dl_elbow_svg <- downloadHandler(
      filename = function() paste0("CoTRA_PCA_elbow_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".svg"),
      content = function(file) save_plot(elbow_plot_gg, file, "svg")
    )
    
    output$dl_pc_variance <- downloadHandler(
      filename = function() paste0("CoTRA_PCA_variance_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
      content = function(file) {
        df <- pc_variance_df()
        req(df)
        utils::write.csv(df, file, row.names = FALSE)
      }
    )
    
    output$dl_pca_scatter_pdf <- downloadHandler(
      filename = function() paste0("CoTRA_PCA_scatter_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"),
      content = function(file) save_plot(scatter_plot_gg, file, "pdf")
    )
    
    output$dl_pca_scatter_svg <- downloadHandler(
      filename = function() paste0("CoTRA_PCA_scatter_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".svg"),
      content = function(file) save_plot(scatter_plot_gg, file, "svg")
    )
    
    output$dl_loadings <- downloadHandler(
      filename = function() paste0("CoTRA_PC_loadings_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
      content = function(file) {
        df <- loading_table_df()
        req(df)
        utils::write.csv(df, file, row.names = FALSE)
      }
    )
    
    output$dl_heatmap_pdf <- downloadHandler(
      filename = function() paste0("CoTRA_PCA_loading_heatmap_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"),
      content = function(file) save_plot(heatmap_plot_gg, file, "pdf", width = 8, height = 10)
    )
    
    output$dl_heatmap_svg <- downloadHandler(
      filename = function() paste0("CoTRA_PCA_loading_heatmap_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".svg"),
      content = function(file) save_plot(heatmap_plot_gg, file, "svg", width = 8, height = 10)
    )
    
    output$download_pca_session <- downloadHandler(
      filename = function() {
        paste0("CoTRA_scRNA_PCA_session_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds")
      },
      content = function(file) {
        seu_out <- pca_seu()
        req(seu_out)
        
        out_dir <- file.path("outputs", "scRNA", "sessioninfo")
        if (!dir.exists(out_dir)) {
          dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
        }
        
        out_file <- file.path(
          out_dir,
          paste0("CoTRA_scRNA_PCA_session_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds")
        )
        
        project_obj <- build_pca_project(seu_out, out_file)
        
        saveRDS(project_obj, out_file)
        file.copy(out_file, file, overwrite = TRUE)
        last_session_path(out_file)
      }
    )
    
    return(reactive({
      if (!is.null(pca_seu())) {
        pca_seu()
      } else {
        get_seu()
      }
    }))
  })
}

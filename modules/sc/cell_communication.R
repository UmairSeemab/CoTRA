# ==========================================================
# modules/sc/cell_communication.R
# CoTRA scRNA-seq Cell Cell Communication Module
#
# Backend:
# - CellChat
#
# Return:
# - seurat
# - comm_ready
# - method
# - results
# - interaction_table
# - session_path
# ==========================================================

suppressPackageStartupMessages({
  library(shiny)
  library(Seurat)
  library(ggplot2)
  library(DT)
})

mod_sc_comm_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    h3("Cell cell communication analysis"),
    
    div(
      class = "alert alert-warning",
      strong("Important: "),
      "Cell cell communication analysis infers ligand receptor activity from expression data. It does not prove physical interaction. Interpret results with cell type annotation, cluster size, sample design, and expression thresholds."
    ),
    
    sidebarLayout(
      sidebarPanel(
        selectInput(
          ns("comm_method"),
          "Communication method",
          choices = c("CellChat" = "cellchat"),
          selected = "cellchat"
        ),
        
        selectInput(
          ns("species"),
          "Species",
          choices = c(
            "Human" = "human",
            "Mouse" = "mouse"
          ),
          selected = "human"
        ),
        
        uiOutput(ns("cluster_col_ui")),
        
        numericInput(
          ns("min_cells"),
          "Minimum cells per cluster",
          value = 10,
          min = 3,
          max = 100,
          step = 1
        ),
        
        numericInput(
          ns("top_n"),
          "Top interactions to show",
          value = 50,
          min = 10,
          max = 500,
          step = 10
        ),
        
        checkboxInput(
          ns("interactive_plots"),
          "Use interactive Plotly plots when possible",
          value = FALSE
        ),
        
        actionButton(
          ns("run_comm"),
          "Run CellChat analysis",
          class = "btn-primary"
        ),
        
        br(), br(),
        downloadButton(ns("download_interactions"), "Download interactions CSV"),
        br(), br(),
        downloadButton(ns("download_summary"), "Download summary CSV"),
        br(), br(),
        downloadButton(ns("download_session"), "Download communication session RDS"),
        br(), br(),
        downloadButton(ns("download_all_zip"), "Download all communication outputs ZIP")
      ),
      
      mainPanel(
        tabsetPanel(
          tabPanel(
            "Summary",
            br(),
            htmlOutput(ns("method_summary")),
            br(),
            DTOutput(ns("summary_table"))
          ),
          
          tabPanel(
            "Heatmap",
            br(),
            uiOutput(ns("heatmap_ui")),
            br(),
            downloadButton(ns("download_heatmap_pdf"), "PDF"),
            downloadButton(ns("download_heatmap_svg"), "SVG")
          ),
          
          tabPanel(
            "Bubble plot",
            br(),
            uiOutput(ns("bubble_ui")),
            br(),
            downloadButton(ns("download_bubble_pdf"), "PDF"),
            downloadButton(ns("download_bubble_svg"), "SVG")
          ),
          
          tabPanel(
            "Network",
            br(),
            uiOutput(ns("network_ui")),
            br(),
            downloadButton(ns("download_network_pdf"), "PDF"),
            downloadButton(ns("download_network_svg"), "SVG")
          ),
          
          tabPanel(
            "Interactions",
            br(),
            DTOutput(ns("interaction_table"))
          ),
          
          tabPanel(
            "Help",
            br(),
            h4("When to use this module"),
            p("Use this module after clustering or annotation when you want to infer likely ligand receptor communication between cell groups."),
            
            h4("Recommended workflow"),
            tags$ul(
              tags$li("Run quality control, normalization, PCA, UMAP, clustering, and markers first."),
              tags$li("Use meaningful cluster labels when available, such as CoTRA_cluster_label."),
              tags$li("Avoid using communication results from very small clusters."),
              tags$li("Validate important ligand receptor pairs with expression plots and known biology.")
            ),
            
            h4("Backend"),
            tags$ul(
              tags$li("CellChat is used as the stable communication backend in this CoTRA module.")
            ),
            
            h4("Interpretation"),
            tags$table(
              class = "table table-bordered table-striped",
              tags$thead(
                tags$tr(
                  tags$th("Output"),
                  tags$th("Meaning"),
                  tags$th("Common mistake")
                )
              ),
              tags$tbody(
                tags$tr(
                  tags$td("Heatmap"),
                  tags$td("Summarizes communication strength between groups."),
                  tags$td("Ignoring group size and expression dropout.")
                ),
                tags$tr(
                  tags$td("Bubble plot"),
                  tags$td("Shows ligand receptor interactions between source and target groups."),
                  tags$td("Treating predicted interaction as direct experimental evidence.")
                ),
                tags$tr(
                  tags$td("Network"),
                  tags$td("Shows global communication links between groups."),
                  tags$td("Overinterpreting weak links from small clusters.")
                ),
                tags$tr(
                  tags$td("Interaction table"),
                  tags$td("Ranked ligand receptor results."),
                  tags$td("Ignoring whether ligand and receptor are expressed in the correct groups.")
                )
              )
            )
          )
        )
      )
    )
  )
}

mod_sc_comm_server <- function(id, seurat_r, sc_state = NULL) {
  moduleServer(id, function(input, output, session) {
    
    ns <- session$ns
    
    `%||%` <- function(a, b) {
      if (!is.null(a)) a else b
    }
    
    get_obj <- reactive({
      obj <- seurat_r()
      validate(need(inherits(obj, "Seurat"), "Input must be a Seurat object."))
      
      if (!"CoTRA_cluster_label" %in% colnames(obj@meta.data) &&
          !"seurat_clusters" %in% colnames(obj@meta.data)) {
        obj$CoTRA_active_ident <- as.character(Idents(obj))
      }
      
      obj
    })
    
    get_assay_data_safe <- function(obj, assay = NULL, layer = "data", slot = "data") {
      assay <- assay %||% DefaultAssay(obj)
      
      out <- tryCatch({
        GetAssayData(obj, assay = assay, layer = layer)
      }, error = function(e) NULL)
      
      if (!is.null(out) && nrow(out) > 0 && ncol(out) > 0) {
        return(out)
      }
      
      out <- tryCatch({
        suppressWarnings(GetAssayData(obj, assay = assay, slot = slot))
      }, error = function(e) NULL)
      
      if (is.null(out) || nrow(out) == 0 || ncol(out) == 0) {
        stop("Could not access normalized assay data.")
      }
      
      out
    }
    
    get_cluster_col <- function(obj) {
      if ("CoTRA_cluster_label" %in% colnames(obj@meta.data)) {
        return("CoTRA_cluster_label")
      }
      
      if ("seurat_clusters" %in% colnames(obj@meta.data)) {
        return("seurat_clusters")
      }
      
      if ("CoTRA_active_ident" %in% colnames(obj@meta.data)) {
        return("CoTRA_active_ident")
      }
      
      NULL
    }
    
    output$cluster_col_ui <- renderUI({
      obj <- get_obj()
      cols <- colnames(obj@meta.data)
      selected <- get_cluster_col(obj)
      
      selectInput(
        ns("cluster_col"),
        "Cluster or cell type column",
        choices = cols,
        selected = selected
      )
    })
    
    comm_res <- reactiveVal(NULL)
    
    check_inputs <- function(obj) {
      validate(need(input$comm_method == "cellchat", "CellChat is the active backend."))
      validate(need(requireNamespace("CellChat", quietly = TRUE), "Install CellChat first."))
      
      if (isTRUE(input$interactive_plots)) {
        validate(need(requireNamespace("plotly", quietly = TRUE), "Install plotly first or disable interactive plots."))
      }
      
      validate(need(!is.null(input$cluster_col), "Select a cluster or cell type column."))
      validate(need(input$cluster_col %in% colnames(obj@meta.data), "Selected cluster column was not found."))
      
      mat_ok <- tryCatch({
        mat <- get_assay_data_safe(obj)
        nrow(mat) > 0 && ncol(mat) > 0
      }, error = function(e) FALSE)
      
      validate(need(mat_ok, "No normalized data found. Run NormalizeData or SCTransform first."))
      
      groups <- as.character(obj@meta.data[[input$cluster_col]])
      groups <- groups[!is.na(groups)]
      
      validate(need(length(unique(groups)) >= 2, "CellChat needs at least two cell groups."))
      
      group_counts <- table(groups)
      valid_groups <- names(group_counts[group_counts >= input$min_cells])
      
      validate(need(length(valid_groups) >= 2, "At least two groups must pass the minimum cell filter."))
      
      TRUE
    }
    
    subset_small_groups <- function(obj) {
      groups <- as.character(obj@meta.data[[input$cluster_col]])
      names(groups) <- rownames(obj@meta.data)
      
      group_counts <- table(groups)
      keep_groups <- names(group_counts[group_counts >= input$min_cells])
      keep_cells <- names(groups)[groups %in% keep_groups]
      
      subset(obj, cells = keep_cells)
    }
    
    run_cellchat <- function(obj) {
      mat <- get_assay_data_safe(obj)
      meta <- obj@meta.data
      
      safe_group_col <- "CoTRA_CellChat_Group"
      
      groups <- as.character(meta[[input$cluster_col]])
      names(groups) <- rownames(meta)
      groups[is.na(groups) | groups == ""] <- "Unknown"
      groups <- paste0("Cluster_", groups)
      groups <- make.names(groups)
      names(groups) <- rownames(meta)
      
      meta[[safe_group_col]] <- factor(groups)
      
      cellchat <- CellChat::createCellChat(
        object = mat,
        meta = meta,
        group.by = safe_group_col
      )
      
      if (input$species == "human") {
        cellchat@DB <- CellChat::CellChatDB.human
      } else {
        cellchat@DB <- CellChat::CellChatDB.mouse
      }
      
      cellchat <- CellChat::subsetData(cellchat)
      cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
      cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)
      cellchat <- CellChat::computeCommunProb(cellchat)
      cellchat <- CellChat::filterCommunication(cellchat, min.cells = input$min_cells)
      cellchat <- CellChat::computeCommunProbPathway(cellchat)
      cellchat <- CellChat::aggregateNet(cellchat)
      
      df <- tryCatch({
        CellChat::subsetCommunication(cellchat)
      }, error = function(e) data.frame())
      
      if (!is.null(df) && nrow(df) > 0 && "prob" %in% colnames(df)) {
        df <- df[order(df$prob, decreasing = TRUE), , drop = FALSE]
      }
      
      mat_count <- tryCatch(cellchat@net$count, error = function(e) NULL)
      mat_weight <- tryCatch(cellchat@net$weight, error = function(e) NULL)
      
      validate(need(!is.null(mat_weight), "CellChat did not return a communication matrix."))
      
      list(
        object = cellchat,
        interaction_table = head(df, input$top_n),
        full_table = df,
        count_matrix = mat_count,
        weight_matrix = mat_weight
      )
    }
    
    matrix_to_long <- function(mat, value_name = "value") {
      if (is.null(mat)) {
        return(data.frame())
      }
      
      df <- as.data.frame(as.table(mat))
      colnames(df) <- c("source", "target", value_name)
      df
    }
    
    observeEvent(input$run_comm, {
      obj <- get_obj()
      check_inputs(obj)
      
      res <- withProgress(message = "Running CellChat analysis", value = 0, {
        
        incProgress(0.2, detail = "Filtering small groups")
        
        obj2 <- subset_small_groups(obj)
        
        incProgress(0.45, detail = "Running CellChat")
        
        backend <- run_cellchat(obj2)
        
        incProgress(0.75, detail = "Saving communication session")
        
        matrix_long <- matrix_to_long(backend$weight_matrix, "communication_strength")
        
        comm_misc <- list(
          method = "cellchat",
          species = input$species,
          cluster_col = input$cluster_col,
          min_cells = input$min_cells,
          top_n = input$top_n,
          interaction_table = backend$interaction_table,
          full_table = backend$full_table,
          count_matrix = backend$count_matrix,
          weight_matrix = backend$weight_matrix,
          matrix_long = matrix_long,
          created = Sys.time()
        )
        
        obj@misc$CoTRA_cell_communication <- comm_misc
        
        dir.create("outputs/scRNA/sessioninfo", recursive = TRUE, showWarnings = FALSE)
        
        timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
        session_path <- file.path(
          "outputs/scRNA/sessioninfo",
          paste0("CoTRA_scRNA_CellCommunication_session_", timestamp, ".rds")
        )
        
        saveRDS(
          list(
            seurat = obj,
            communication = comm_misc,
            communication_object = backend$object,
            sessionInfo = sessionInfo()
          ),
          session_path
        )
        
        incProgress(1, detail = "Done")
        
        list(
          seurat = obj,
          comm_ready = TRUE,
          method = "cellchat",
          species = input$species,
          cluster_col = input$cluster_col,
          min_cells = input$min_cells,
          results = backend$object,
          interaction_table = backend$interaction_table,
          full_table = backend$full_table,
          count_matrix = backend$count_matrix,
          weight_matrix = backend$weight_matrix,
          matrix_long = matrix_long,
          session_path = session_path
        )
      })
      
      comm_res(res)
    })
    
    get_res <- reactive({
      res <- comm_res()
      validate(need(!is.null(res), "Run CellChat analysis first."))
      res
    })
    
    heatmap_plot_obj <- reactive({
      res <- get_res()
      
      df <- res$matrix_long
      validate(need(nrow(df) > 0, "No communication matrix available for heatmap."))
      
      ggplot(df, aes(source, target, fill = communication_strength)) +
        geom_tile(color = "white") +
        scale_fill_viridis_c() +
        theme_classic(base_size = 13) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(
          title = "CellChat communication heatmap",
          x = "Source",
          y = "Target",
          fill = "Strength"
        )
    })
    
    bubble_plot_obj <- reactive({
      res <- get_res()
      df <- res$interaction_table
      
      validate(need(nrow(df) > 0, "No interaction table available for bubble plot."))
      
      source_col <- intersect(c("source", "sender", "source_celltype", "source_label"), colnames(df))[1]
      target_col <- intersect(c("target", "receiver", "target_celltype", "target_label"), colnames(df))[1]
      ligand_col <- intersect(c("ligand", "ligand.complex", "ligand_complex"), colnames(df))[1]
      receptor_col <- intersect(c("receptor", "receptor.complex", "receptor_complex"), colnames(df))[1]
      score_col <- intersect(c("prob", "communication_score", "lrscore", "score"), colnames(df))[1]
      
      validate(need(!is.na(source_col) && !is.na(target_col), "Source and target columns were not found."))
      
      if (is.na(ligand_col)) df$ligand_plot <- "ligand" else df$ligand_plot <- as.character(df[[ligand_col]])
      if (is.na(receptor_col)) df$receptor_plot <- "receptor" else df$receptor_plot <- as.character(df[[receptor_col]])
      
      if (is.na(score_col)) {
        df$score_plot <- seq_len(nrow(df))
      } else {
        df$score_plot <- suppressWarnings(as.numeric(df[[score_col]]))
      }
      
      df$pair <- paste(df$ligand_plot, df$receptor_plot, sep = " -> ")
      df$source_target <- paste(df[[source_col]], df[[target_col]], sep = " -> ")
      
      ggplot(df, aes(source_target, pair)) +
        geom_point(aes(size = abs(score_plot), color = score_plot), alpha = 0.85) +
        scale_color_viridis_c() +
        theme_classic(base_size = 12) +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 8)
        ) +
        labs(
          title = "Top CellChat interactions",
          x = "Source -> Target",
          y = "Ligand -> Receptor",
          size = "Score",
          color = "Score"
        )
    })
    
    network_plot_obj <- reactive({
      res <- get_res()
      
      df <- res$matrix_long
      validate(need(nrow(df) > 0, "No communication matrix available for network plot."))
      
      df <- df[df$communication_strength > 0, , drop = FALSE]
      validate(need(nrow(df) > 0, "No positive communication links found."))
      
      ggplot(df, aes(source, target)) +
        geom_point(aes(size = communication_strength, color = communication_strength), alpha = 0.85) +
        scale_color_viridis_c() +
        theme_classic(base_size = 12) +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 8)
        ) +
        labs(
          title = "CellChat communication network summary",
          x = "Source",
          y = "Target",
          size = "Strength",
          color = "Strength"
        )
    })
    
    output$heatmap_ui <- renderUI({
      if (isTRUE(input$interactive_plots)) {
        plotly::plotlyOutput(ns("heatmap_plotly"), height = "650px")
      } else {
        plotOutput(ns("heatmap"), height = "650px")
      }
    })
    
    output$bubble_ui <- renderUI({
      if (isTRUE(input$interactive_plots)) {
        plotly::plotlyOutput(ns("bubble_plotly"), height = "650px")
      } else {
        plotOutput(ns("bubble"), height = "650px")
      }
    })
    
    output$network_ui <- renderUI({
      if (isTRUE(input$interactive_plots)) {
        plotly::plotlyOutput(ns("network_plotly"), height = "650px")
      } else {
        plotOutput(ns("network"), height = "650px")
      }
    })
    
    output$heatmap <- renderPlot({ heatmap_plot_obj() })
    output$bubble <- renderPlot({ bubble_plot_obj() })
    output$network <- renderPlot({ network_plot_obj() })
    
    if (requireNamespace("plotly", quietly = TRUE)) {
      output$heatmap_plotly <- plotly::renderPlotly({
        plotly::ggplotly(heatmap_plot_obj())
      })
      
      output$bubble_plotly <- plotly::renderPlotly({
        plotly::ggplotly(bubble_plot_obj())
      })
      
      output$network_plotly <- plotly::renderPlotly({
        plotly::ggplotly(network_plot_obj())
      })
    }
    
    output$method_summary <- renderUI({
      res <- get_res()
      
      HTML(paste0(
        "<div class='alert alert-info'>",
        "<b>Method:</b> CellChat<br>",
        "<b>Species:</b> ", res$species, "<br>",
        "<b>Cluster column:</b> ", res$cluster_col, "<br>",
        "<b>Minimum cells per group:</b> ", res$min_cells, "<br>",
        "<b>Interactions shown:</b> ", nrow(res$interaction_table), "<br>",
        "<b>Total interactions:</b> ", nrow(res$full_table), "<br>",
        "<b>Session path:</b> ", res$session_path,
        "</div>"
      ))
    })
    
    output$summary_table <- renderDT({
      res <- get_res()
      
      summary_df <- data.frame(
        metric = c(
          "Method",
          "Species",
          "Cluster column",
          "Minimum cells per group",
          "Interactions shown",
          "Total interactions",
          "Session path"
        ),
        value = c(
          "CellChat",
          res$species,
          res$cluster_col,
          res$min_cells,
          nrow(res$interaction_table),
          nrow(res$full_table),
          res$session_path
        ),
        stringsAsFactors = FALSE
      )
      
      datatable(summary_df, options = list(dom = "t", paging = FALSE), rownames = FALSE)
    })
    
    output$interaction_table <- renderDT({
      res <- get_res()
      
      datatable(
        res$interaction_table,
        options = list(pageLength = 20, scrollX = TRUE),
        rownames = FALSE
      )
    })
    
    save_plot <- function(plot_obj, file, width = 10, height = 8) {
      ext <- tools::file_ext(file)
      
      if (tolower(ext) == "svg") {
        validate(need(requireNamespace("svglite", quietly = TRUE), "Install svglite to export SVG files."))
      }
      
      ggplot2::ggsave(
        filename = file,
        plot = plot_obj,
        width = width,
        height = height,
        units = "in",
        limitsize = FALSE
      )
    }
    
    output$download_interactions <- downloadHandler(
      filename = function() paste0("CoTRA_cellchat_interactions_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
      content = function(file) {
        write.csv(get_res()$full_table, file, row.names = FALSE)
      }
    )
    
    output$download_summary <- downloadHandler(
      filename = function() paste0("CoTRA_cellchat_summary_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
      content = function(file) {
        res <- get_res()
        summary_df <- data.frame(
          metric = c("method", "species", "cluster_col", "min_cells", "interactions_shown", "total_interactions", "session_path"),
          value = c("cellchat", res$species, res$cluster_col, res$min_cells, nrow(res$interaction_table), nrow(res$full_table), res$session_path)
        )
        write.csv(summary_df, file, row.names = FALSE)
      }
    )
    
    output$download_session <- downloadHandler(
      filename = function() basename(get_res()$session_path),
      content = function(file) {
        file.copy(get_res()$session_path, file, overwrite = TRUE)
      }
    )
    
    output$download_heatmap_pdf <- downloadHandler(
      filename = function() paste0("CoTRA_cellchat_heatmap_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"),
      content = function(file) save_plot(heatmap_plot_obj(), file)
    )
    
    output$download_heatmap_svg <- downloadHandler(
      filename = function() paste0("CoTRA_cellchat_heatmap_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".svg"),
      content = function(file) save_plot(heatmap_plot_obj(), file)
    )
    
    output$download_bubble_pdf <- downloadHandler(
      filename = function() paste0("CoTRA_cellchat_bubble_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"),
      content = function(file) save_plot(bubble_plot_obj(), file, width = 12, height = 9)
    )
    
    output$download_bubble_svg <- downloadHandler(
      filename = function() paste0("CoTRA_cellchat_bubble_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".svg"),
      content = function(file) save_plot(bubble_plot_obj(), file, width = 12, height = 9)
    )
    
    output$download_network_pdf <- downloadHandler(
      filename = function() paste0("CoTRA_cellchat_network_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"),
      content = function(file) save_plot(network_plot_obj(), file)
    )
    
    output$download_network_svg <- downloadHandler(
      filename = function() paste0("CoTRA_cellchat_network_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".svg"),
      content = function(file) save_plot(network_plot_obj(), file)
    )
    
    output$download_all_zip <- downloadHandler(
      filename = function() paste0("CoTRA_cellchat_outputs_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".zip"),
      content = function(file) {
        res <- get_res()
        
        tmpdir <- tempfile("CoTRA_cellchat_outputs_")
        dir.create(tmpdir, recursive = TRUE, showWarnings = FALSE)
        
        write.csv(res$full_table, file.path(tmpdir, "CoTRA_cellchat_interactions.csv"), row.names = FALSE)
        write.csv(res$matrix_long, file.path(tmpdir, "CoTRA_cellchat_matrix_long.csv"), row.names = FALSE)
        
        file.copy(res$session_path, file.path(tmpdir, basename(res$session_path)), overwrite = TRUE)
        
        save_plot(heatmap_plot_obj(), file.path(tmpdir, "CoTRA_cellchat_heatmap.pdf"))
        save_plot(bubble_plot_obj(), file.path(tmpdir, "CoTRA_cellchat_bubble.pdf"), width = 12, height = 9)
        save_plot(network_plot_obj(), file.path(tmpdir, "CoTRA_cellchat_network.pdf"))
        
        if (requireNamespace("svglite", quietly = TRUE)) {
          save_plot(heatmap_plot_obj(), file.path(tmpdir, "CoTRA_cellchat_heatmap.svg"))
          save_plot(bubble_plot_obj(), file.path(tmpdir, "CoTRA_cellchat_bubble.svg"), width = 12, height = 9)
          save_plot(network_plot_obj(), file.path(tmpdir, "CoTRA_cellchat_network.svg"))
        }
        
        oldwd <- getwd()
        on.exit(setwd(oldwd), add = TRUE)
        setwd(tmpdir)
        
        utils::zip(zipfile = file, files = list.files(tmpdir, recursive = FALSE))
      }
    )
    
    return(list(
      seurat = reactive({
        res <- get_res()
        res$seurat
      }),
      
      comm_ready = reactive({
        res <- get_res()
        isTRUE(res$comm_ready)
      }),
      
      method = reactive({
        "cellchat"
      }),
      
      results = reactive({
        res <- get_res()
        res$results
      }),
      
      interaction_table = reactive({
        res <- get_res()
        res$interaction_table
      }),
      
      session_path = reactive({
        res <- get_res()
        res$session_path
      })
    ))
  })
}
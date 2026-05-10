# ==========================================================
# modules/sc/differential_abundance.R
# CoTRA scRNA-seq Differential Abundance Module
#
# Methods:
# - propeller-like limma
# - edgeR
# - Fisher exact test
# - Chi-square test
#
# Outputs:
# - DA summary
# - Count table
# - Proportion table
# - DA result table
# - Stacked proportion plot
# - Boxplot per cell type
# - Volcano plot
# - LogFC heatmap
# - UMAP plot
# - PDF and SVG plot downloads
# - HTML interactive plot download
# - ZIP export
# - Session RDS save
#
# Return:
# - seurat
# - da_ready
# - da_counts
# - da_proportions
# - da_results
# - method
# - session_path
# ==========================================================

suppressPackageStartupMessages({
  library(shiny)
  library(Seurat)
  library(ggplot2)
  library(DT)
  library(dplyr)
  library(tidyr)
  library(scales)
})

mod_sc_differential_abundance_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    h3("Differential abundance analysis"),
    
    div(
      class = "alert alert-warning",
      strong("Important: "),
      "Differential abundance tests whether cell type or cluster proportions differ between biological conditions. Use biological samples as replicates, not individual cells."
    ),
    
    sidebarLayout(
      sidebarPanel(
        selectInput(
          ns("method"),
          "Differential abundance method",
          choices = c(
            "Auto" = "auto",
            "propeller-like limma" = "propeller",
            "edgeR" = "edger",
            "Fisher exact test" = "fisher",
            "Chi-square test" = "chisq"
          ),
          selected = "auto"
        ),
        
        uiOutput(ns("group_col_ui")),
        uiOutput(ns("sample_col_ui")),
        uiOutput(ns("celltype_col_ui")),
        uiOutput(ns("contrast_ui")),
        
        numericInput(
          ns("min_cells"),
          "Minimum cells per cell type",
          value = 20,
          min = 1,
          step = 1
        ),
        
        numericInput(
          ns("fdr_cutoff"),
          "FDR cutoff",
          value = 0.05,
          min = 0,
          max = 1,
          step = 0.01
        ),
        
        checkboxInput(
          ns("interactive_plots"),
          "Use interactive Plotly plots when possible",
          value = FALSE
        ),
        
        actionButton(
          ns("run_da"),
          "Run differential abundance",
          class = "btn-primary"
        ),
        
        br(), br(),
        downloadButton(ns("download_da_csv"), "Download DA results CSV"),
        br(), br(),
        downloadButton(ns("download_counts_csv"), "Download counts CSV"),
        br(), br(),
        downloadButton(ns("download_props_csv"), "Download proportions CSV"),
        br(), br(),
        downloadButton(ns("download_session"), "Download DA session RDS"),
        br(), br(),
        downloadButton(ns("download_all_zip"), "Download all DA outputs ZIP")
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
            "DA results",
            br(),
            DTOutput(ns("da_table"))
          ),
          
          tabPanel(
            "Counts",
            br(),
            DTOutput(ns("count_table"))
          ),
          
          tabPanel(
            "Proportions",
            br(),
            DTOutput(ns("prop_table"))
          ),
          
          tabPanel(
            "Stacked proportions",
            br(),
            uiOutput(ns("stacked_plot_ui")),
            br(),
            downloadButton(ns("download_stacked_pdf"), "PDF"),
            downloadButton(ns("download_stacked_svg"), "SVG"),
            downloadButton(ns("download_stacked_html"), "HTML")
          ),
          
          tabPanel(
            "Cell type boxplot",
            br(),
            uiOutput(ns("boxplot_celltype_ui")),
            br(),
            uiOutput(ns("boxplot_ui")),
            br(),
            downloadButton(ns("download_boxplot_pdf"), "PDF"),
            downloadButton(ns("download_boxplot_svg"), "SVG"),
            downloadButton(ns("download_boxplot_html"), "HTML")
          ),
          
          tabPanel(
            "Volcano plot",
            br(),
            uiOutput(ns("volcano_plot_ui")),
            br(),
            downloadButton(ns("download_volcano_pdf"), "PDF"),
            downloadButton(ns("download_volcano_svg"), "SVG"),
            downloadButton(ns("download_volcano_html"), "HTML")
          ),
          
          tabPanel(
            "LogFC heatmap",
            br(),
            uiOutput(ns("heatmap_plot_ui")),
            br(),
            downloadButton(ns("download_heatmap_pdf"), "PDF"),
            downloadButton(ns("download_heatmap_svg"), "SVG"),
            downloadButton(ns("download_heatmap_html"), "HTML")
          ),
          
          tabPanel(
            "UMAP",
            br(),
            uiOutput(ns("umap_plot_ui")),
            br(),
            downloadButton(ns("download_umap_pdf"), "PDF"),
            downloadButton(ns("download_umap_svg"), "SVG"),
            downloadButton(ns("download_umap_html"), "HTML")
          ),
          
          tabPanel(
            "Help",
            br(),
            h4("When to use this module"),
            p("Use this module after annotation when you want to test whether cell type abundance differs between conditions."),
            
            h4("Required columns"),
            tags$ul(
              tags$li("Condition column: disease, control, treatment, genotype, time point."),
              tags$li("Sample column: biological sample, donor, replicate, orig.ident."),
              tags$li("Cell type column: CoTRA_celltype, annotation, cell_type, seurat_clusters.")
            ),
            
            h4("Methods"),
            tags$ul(
              tags$li("propeller-like limma uses sample-level cell type proportions with asin-square-root transformation."),
              tags$li("edgeR models cell type counts across samples."),
              tags$li("Fisher and chi-square are fallback tests on total cell counts.")
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
                  tags$td("Stacked proportions"),
                  tags$td("Cell composition per sample."),
                  tags$td("Using this plot without checking replicate variability.")
                ),
                tags$tr(
                  tags$td("Boxplot"),
                  tags$td("Sample-level proportion distribution per group."),
                  tags$td("Treating individual cells as replicates.")
                ),
                tags$tr(
                  tags$td("Volcano plot"),
                  tags$td("Effect size and significance of abundance changes."),
                  tags$td("Focusing only on FDR without checking effect size.")
                ),
                tags$tr(
                  tags$td("LogFC heatmap"),
                  tags$td("Direction and size of abundance changes."),
                  tags$td("Overinterpreting rare cell types with low cell counts.")
                )
              )
            )
          )
        )
      )
    )
  )
}

mod_sc_differential_abundance_server <- function(
    id,
    seurat_r,
    sc_state = NULL,
    output_dir = "outputs/scRNA/differential_abundance"
) {
  moduleServer(id, function(input, output, session) {
    
    ns <- session$ns
    da_res <- reactiveVal(NULL)
    
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(output_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(output_dir, "sessioninfo"), recursive = TRUE, showWarnings = FALSE)
    
    get_obj <- reactive({
      obj <- seurat_r()
      validate(need(inherits(obj, "Seurat"), "Input must be a Seurat object."))
      obj
    })
    
    get_meta <- reactive({
      obj <- get_obj()
      md <- obj@meta.data
      md$cell_id <- rownames(md)
      md
    })
    
    guess_sc_da_column <- function(cols, patterns) {
      for (p in patterns) {
        hit <- cols[grepl(p, cols, ignore.case = TRUE)]
        if (length(hit) > 0) return(hit[1])
      }
      cols[1]
    }
    
    output$group_col_ui <- renderUI({
      md <- get_meta()
      cols <- colnames(md)
      selected <- guess_sc_da_column(cols, c("condition", "group", "treatment", "status", "genotype"))
      
      selectInput(
        ns("group_col"),
        "Condition column",
        choices = cols,
        selected = selected
      )
    })
    
    output$sample_col_ui <- renderUI({
      md <- get_meta()
      cols <- colnames(md)
      selected <- guess_sc_da_column(cols, c("sample", "sample_id", "orig.ident", "donor", "replicate"))
      
      selectInput(
        ns("sample_col"),
        "Sample / replicate column",
        choices = cols,
        selected = selected
      )
    })
    
    output$celltype_col_ui <- renderUI({
      md <- get_meta()
      cols <- colnames(md)
      selected <- guess_sc_da_column(
        cols,
        c("CoTRA_celltype", "celltype", "cell_type", "annotation", "predicted", "seurat_clusters", "cluster")
      )
      
      selectInput(
        ns("celltype_col"),
        "Cell type / cluster column",
        choices = cols,
        selected = selected
      )
    })
    
    output$contrast_ui <- renderUI({
      req(input$group_col)
      md <- get_meta()
      
      groups <- sort(unique(as.character(md[[input$group_col]])))
      groups <- groups[!is.na(groups) & groups != ""]
      
      if (length(groups) < 2) {
        return(helpText("At least two groups are required."))
      }
      
      tagList(
        selectInput(
          ns("group_1"),
          "Contrast group 1",
          choices = groups,
          selected = groups[1]
        ),
        selectInput(
          ns("group_2"),
          "Contrast group 2",
          choices = groups,
          selected = groups[2]
        )
      )
    })
    
    output$boxplot_celltype_ui <- renderUI({
      res <- get_res()
      celltypes <- sort(unique(as.character(res$long$.celltype)))
      
      selectInput(
        ns("selected_celltype"),
        "Cell type for boxplot",
        choices = celltypes,
        selected = celltypes[1]
      )
    })
    
    choose_sc_da_method <- function(method) {
      if (method != "auto") return(method)
      if (requireNamespace("limma", quietly = TRUE)) return("propeller")
      if (requireNamespace("edgeR", quietly = TRUE)) return("edger")
      "fisher"
    }
    
    observeEvent(input$run_da, {
      obj <- get_obj()
      md <- get_meta()
      
      validate(
        need(input$group_col %in% colnames(md), "Selected condition column was not found."),
        need(input$sample_col %in% colnames(md), "Selected sample column was not found."),
        need(input$celltype_col %in% colnames(md), "Selected cell type column was not found."),
        need(input$group_1 != input$group_2, "Contrast groups must be different.")
      )
      
      res <- withProgress(message = "Running differential abundance analysis", value = 0, {
        
        incProgress(0.2, detail = "Preparing sample-level count table")
        
        md2 <- md %>%
          mutate(
            .group = as.character(.data[[input$group_col]]),
            .sample = as.character(.data[[input$sample_col]]),
            .celltype = as.character(.data[[input$celltype_col]])
          ) %>%
          filter(
            !is.na(.group),
            !is.na(.sample),
            !is.na(.celltype),
            .group %in% c(input$group_1, input$group_2),
            .group != "",
            .sample != "",
            .celltype != ""
          )
        
        validate(
          need(nrow(md2) > 0, "No cells found for this contrast."),
          need(length(unique(md2$.group)) == 2, "Selected data must contain two groups."),
          need(length(unique(md2$.sample)) >= 2, "At least two samples are required.")
        )
        
        valid_celltypes <- md2 %>%
          count(.celltype, name = "total_cells") %>%
          filter(total_cells >= input$min_cells) %>%
          pull(.celltype)
        
        md2 <- md2 %>% filter(.celltype %in% valid_celltypes)
        
        validate(
          need(nrow(md2) > 0, "No cell types passed the minimum cell filter.")
        )
        
        counts_long <- md2 %>%
          count(.sample, .group, .celltype, name = "n_cells") %>%
          complete(.sample, .group, .celltype, fill = list(n_cells = 0)) %>%
          group_by(.sample) %>%
          mutate(
            total_sample_cells = sum(n_cells),
            proportion = ifelse(total_sample_cells > 0, n_cells / total_sample_cells, 0)
          ) %>%
          ungroup()
        
        counts_wide <- counts_long %>%
          select(.sample, .group, .celltype, n_cells) %>%
          pivot_wider(
            names_from = .celltype,
            values_from = n_cells,
            values_fill = 0
          )
        
        props_wide <- counts_long %>%
          select(.sample, .group, .celltype, proportion) %>%
          pivot_wider(
            names_from = .celltype,
            values_from = proportion,
            values_fill = 0
          )
        
        incProgress(0.45, detail = "Running selected method")
        
        method_used <- choose_sc_da_method(input$method)
        
        results <- run_sc_da_method(
          counts = counts_long,
          method = method_used,
          group_1 = input$group_1,
          group_2 = input$group_2
        ) %>%
          mutate(
            significant = ifelse(FDR <= input$fdr_cutoff, "yes", "no"),
            method = method_used
          ) %>%
          arrange(FDR)
        
        incProgress(0.75, detail = "Saving DA session")
        
        timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
        session_path <- file.path(
          output_dir,
          "sessioninfo",
          paste0("CoTRA_scRNA_DifferentialAbundance_session_", timestamp, ".rds")
        )
        
        out <- list(
          seurat = obj,
          da_ready = TRUE,
          method = method_used,
          group_col = input$group_col,
          sample_col = input$sample_col,
          celltype_col = input$celltype_col,
          group_1 = input$group_1,
          group_2 = input$group_2,
          metadata = md2,
          long = counts_long,
          counts = counts_wide,
          proportions = props_wide,
          results = results,
          session_path = session_path,
          created = Sys.time()
        )
        
        saveRDS(
          list(
            seurat_metadata = md,
            differential_abundance = out,
            sessionInfo = sessionInfo()
          ),
          session_path
        )
        
        obj@misc$CoTRA_differential_abundance <- out
        
        if (!is.null(sc_state) && is.environment(sc_state)) {
          sc_state$da <- out
        }
        
        incProgress(1, detail = "Done")
        
        out
      })
      
      da_res(res)
    })
    
    get_res <- reactive({
      res <- da_res()
      validate(need(!is.null(res), "Run differential abundance analysis first."))
      res
    })
    
    output$method_summary <- renderUI({
      res <- get_res()
      
      HTML(paste0(
        "<div class='alert alert-info'>",
        "<b>Method:</b> ", res$method, "<br>",
        "<b>Condition column:</b> ", res$group_col, "<br>",
        "<b>Sample column:</b> ", res$sample_col, "<br>",
        "<b>Cell type column:</b> ", res$celltype_col, "<br>",
        "<b>Contrast:</b> ", res$group_1, " vs ", res$group_2, "<br>",
        "<b>Cell types tested:</b> ", nrow(res$results), "<br>",
        "<b>Significant cell types:</b> ", sum(res$results$FDR <= input$fdr_cutoff, na.rm = TRUE), "<br>",
        "<b>Session path:</b> ", res$session_path,
        "</div>"
      ))
    })
    
    output$summary_table <- renderDT({
      res <- get_res()
      
      summary_df <- data.frame(
        metric = c(
          "Method",
          "Condition column",
          "Sample column",
          "Cell type column",
          "Contrast group 1",
          "Contrast group 2",
          "Cell types tested",
          "Significant cell types",
          "Session path"
        ),
        value = c(
          res$method,
          res$group_col,
          res$sample_col,
          res$celltype_col,
          res$group_1,
          res$group_2,
          nrow(res$results),
          sum(res$results$FDR <= input$fdr_cutoff, na.rm = TRUE),
          res$session_path
        ),
        stringsAsFactors = FALSE
      )
      
      datatable(summary_df, options = list(dom = "t", paging = FALSE), rownames = FALSE)
    })
    
    output$da_table <- renderDT({
      res <- get_res()
      datatable(res$results, options = list(pageLength = 20, scrollX = TRUE), rownames = FALSE)
    })
    
    output$count_table <- renderDT({
      res <- get_res()
      datatable(res$counts, options = list(pageLength = 20, scrollX = TRUE), rownames = FALSE)
    })
    
    output$prop_table <- renderDT({
      res <- get_res()
      datatable(res$proportions, options = list(pageLength = 20, scrollX = TRUE), rownames = FALSE)
    })
    
    stacked_plot_obj <- reactive({
      res <- get_res()
      plot_sc_da_stacked(res$long)
    })
    
    boxplot_obj <- reactive({
      res <- get_res()
      req(input$selected_celltype)
      plot_sc_da_boxplot(res$long, input$selected_celltype)
    })
    
    volcano_plot_obj <- reactive({
      res <- get_res()
      plot_sc_da_volcano(res$results, input$fdr_cutoff)
    })
    
    heatmap_plot_obj <- reactive({
      res <- get_res()
      plot_sc_da_heatmap(res$results)
    })
    
    umap_plot_obj <- reactive({
      res <- get_res()
      plot_sc_da_umap(res$seurat, res$celltype_col)
    })
    
    output$stacked_plot_ui <- renderUI({
      if (isTRUE(input$interactive_plots)) {
        plotly::plotlyOutput(ns("stacked_plotly"), height = "650px")
      } else {
        plotOutput(ns("stacked_plot"), height = "650px")
      }
    })
    
    output$boxplot_ui <- renderUI({
      if (isTRUE(input$interactive_plots)) {
        plotly::plotlyOutput(ns("boxplot_plotly"), height = "650px")
      } else {
        plotOutput(ns("boxplot_plot"), height = "650px")
      }
    })
    
    output$volcano_plot_ui <- renderUI({
      if (isTRUE(input$interactive_plots)) {
        plotly::plotlyOutput(ns("volcano_plotly"), height = "650px")
      } else {
        plotOutput(ns("volcano_plot"), height = "650px")
      }
    })
    
    output$heatmap_plot_ui <- renderUI({
      if (isTRUE(input$interactive_plots)) {
        plotly::plotlyOutput(ns("heatmap_plotly"), height = "650px")
      } else {
        plotOutput(ns("heatmap_plot"), height = "650px")
      }
    })
    
    output$umap_plot_ui <- renderUI({
      if (isTRUE(input$interactive_plots)) {
        plotly::plotlyOutput(ns("umap_plotly"), height = "650px")
      } else {
        plotOutput(ns("umap_plot"), height = "650px")
      }
    })
    
    output$stacked_plot <- renderPlot({ stacked_plot_obj() })
    output$boxplot_plot <- renderPlot({ boxplot_obj() })
    output$volcano_plot <- renderPlot({ volcano_plot_obj() })
    output$heatmap_plot <- renderPlot({ heatmap_plot_obj() })
    output$umap_plot <- renderPlot({ umap_plot_obj() })
    
    if (requireNamespace("plotly", quietly = TRUE)) {
      output$stacked_plotly <- plotly::renderPlotly({
        plotly::ggplotly(stacked_plot_obj(), tooltip = "text")
      })
      
      output$boxplot_plotly <- plotly::renderPlotly({
        plotly::ggplotly(boxplot_obj(), tooltip = "text")
      })
      
      output$volcano_plotly <- plotly::renderPlotly({
        plotly::ggplotly(volcano_plot_obj(), tooltip = "text")
      })
      
      output$heatmap_plotly <- plotly::renderPlotly({
        plotly::ggplotly(heatmap_plot_obj(), tooltip = "text")
      })
      
      output$umap_plotly <- plotly::renderPlotly({
        plotly::ggplotly(umap_plot_obj(), tooltip = "text")
      })
    }
    
    save_plot <- function(plot_obj, file, width = 10, height = 8) {
      ext <- tolower(tools::file_ext(file))
      
      if (ext == "svg") {
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
    
    save_html_plot <- function(plot_obj, file) {
      validate(need(requireNamespace("plotly", quietly = TRUE), "Install plotly to export HTML plots."))
      validate(need(requireNamespace("htmlwidgets", quietly = TRUE), "Install htmlwidgets to export HTML plots."))
      
      p <- plotly::ggplotly(plot_obj, tooltip = "text")
      htmlwidgets::saveWidget(p, file = file, selfcontained = TRUE)
    }
    
    output$download_da_csv <- downloadHandler(
      filename = function() paste0("CoTRA_scDA_results_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
      content = function(file) write.csv(get_res()$results, file, row.names = FALSE)
    )
    
    output$download_counts_csv <- downloadHandler(
      filename = function() paste0("CoTRA_scDA_counts_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
      content = function(file) write.csv(get_res()$counts, file, row.names = FALSE)
    )
    
    output$download_props_csv <- downloadHandler(
      filename = function() paste0("CoTRA_scDA_proportions_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
      content = function(file) write.csv(get_res()$proportions, file, row.names = FALSE)
    )
    
    output$download_session <- downloadHandler(
      filename = function() basename(get_res()$session_path),
      content = function(file) file.copy(get_res()$session_path, file, overwrite = TRUE)
    )
    
    output$download_stacked_pdf <- downloadHandler(
      filename = function() paste0("CoTRA_scDA_stacked_proportions_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"),
      content = function(file) save_plot(stacked_plot_obj(), file)
    )
    
    output$download_stacked_svg <- downloadHandler(
      filename = function() paste0("CoTRA_scDA_stacked_proportions_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".svg"),
      content = function(file) save_plot(stacked_plot_obj(), file)
    )
    
    output$download_stacked_html <- downloadHandler(
      filename = function() paste0("CoTRA_scDA_stacked_proportions_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".html"),
      content = function(file) save_html_plot(stacked_plot_obj(), file)
    )
    
    output$download_boxplot_pdf <- downloadHandler(
      filename = function() paste0("CoTRA_scDA_boxplot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"),
      content = function(file) save_plot(boxplot_obj(), file)
    )
    
    output$download_boxplot_svg <- downloadHandler(
      filename = function() paste0("CoTRA_scDA_boxplot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".svg"),
      content = function(file) save_plot(boxplot_obj(), file)
    )
    
    output$download_boxplot_html <- downloadHandler(
      filename = function() paste0("CoTRA_scDA_boxplot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".html"),
      content = function(file) save_html_plot(boxplot_obj(), file)
    )
    
    output$download_volcano_pdf <- downloadHandler(
      filename = function() paste0("CoTRA_scDA_volcano_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"),
      content = function(file) save_plot(volcano_plot_obj(), file)
    )
    
    output$download_volcano_svg <- downloadHandler(
      filename = function() paste0("CoTRA_scDA_volcano_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".svg"),
      content = function(file) save_plot(volcano_plot_obj(), file)
    )
    
    output$download_volcano_html <- downloadHandler(
      filename = function() paste0("CoTRA_scDA_volcano_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".html"),
      content = function(file) save_html_plot(volcano_plot_obj(), file)
    )
    
    output$download_heatmap_pdf <- downloadHandler(
      filename = function() paste0("CoTRA_scDA_heatmap_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"),
      content = function(file) save_plot(heatmap_plot_obj(), file)
    )
    
    output$download_heatmap_svg <- downloadHandler(
      filename = function() paste0("CoTRA_scDA_heatmap_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".svg"),
      content = function(file) save_plot(heatmap_plot_obj(), file)
    )
    
    output$download_heatmap_html <- downloadHandler(
      filename = function() paste0("CoTRA_scDA_heatmap_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".html"),
      content = function(file) save_html_plot(heatmap_plot_obj(), file)
    )
    
    output$download_umap_pdf <- downloadHandler(
      filename = function() paste0("CoTRA_scDA_umap_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"),
      content = function(file) save_plot(umap_plot_obj(), file)
    )
    
    output$download_umap_svg <- downloadHandler(
      filename = function() paste0("CoTRA_scDA_umap_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".svg"),
      content = function(file) save_plot(umap_plot_obj(), file)
    )
    
    output$download_umap_html <- downloadHandler(
      filename = function() paste0("CoTRA_scDA_umap_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".html"),
      content = function(file) save_html_plot(umap_plot_obj(), file)
    )
    
    output$download_all_zip <- downloadHandler(
      filename = function() paste0("CoTRA_scDA_outputs_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".zip"),
      content = function(file) {
        res <- get_res()
        
        tmpdir <- tempfile("CoTRA_scDA_outputs_")
        dir.create(tmpdir, recursive = TRUE, showWarnings = FALSE)
        
        write.csv(res$results, file.path(tmpdir, "CoTRA_scDA_results.csv"), row.names = FALSE)
        write.csv(res$counts, file.path(tmpdir, "CoTRA_scDA_counts.csv"), row.names = FALSE)
        write.csv(res$proportions, file.path(tmpdir, "CoTRA_scDA_proportions.csv"), row.names = FALSE)
        write.csv(res$long, file.path(tmpdir, "CoTRA_scDA_long_table.csv"), row.names = FALSE)
        
        file.copy(res$session_path, file.path(tmpdir, basename(res$session_path)), overwrite = TRUE)
        
        save_plot(stacked_plot_obj(), file.path(tmpdir, "CoTRA_scDA_stacked_proportions.pdf"))
        save_plot(boxplot_obj(), file.path(tmpdir, "CoTRA_scDA_boxplot.pdf"))
        save_plot(volcano_plot_obj(), file.path(tmpdir, "CoTRA_scDA_volcano.pdf"))
        save_plot(heatmap_plot_obj(), file.path(tmpdir, "CoTRA_scDA_heatmap.pdf"))
        save_plot(umap_plot_obj(), file.path(tmpdir, "CoTRA_scDA_umap.pdf"))
        
        if (requireNamespace("svglite", quietly = TRUE)) {
          save_plot(stacked_plot_obj(), file.path(tmpdir, "CoTRA_scDA_stacked_proportions.svg"))
          save_plot(boxplot_obj(), file.path(tmpdir, "CoTRA_scDA_boxplot.svg"))
          save_plot(volcano_plot_obj(), file.path(tmpdir, "CoTRA_scDA_volcano.svg"))
          save_plot(heatmap_plot_obj(), file.path(tmpdir, "CoTRA_scDA_heatmap.svg"))
          save_plot(umap_plot_obj(), file.path(tmpdir, "CoTRA_scDA_umap.svg"))
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
      
      da_ready = reactive({
        res <- get_res()
        isTRUE(res$da_ready)
      }),
      
      da_counts = reactive({
        res <- get_res()
        res$counts
      }),
      
      da_proportions = reactive({
        res <- get_res()
        res$proportions
      }),
      
      da_results = reactive({
        res <- get_res()
        res$results
      }),
      
      method = reactive({
        res <- get_res()
        res$method
      }),
      
      session_path = reactive({
        res <- get_res()
        res$session_path
      })
    ))
  })
}

# ==========================================================
# Statistical helper functions
# ==========================================================

run_sc_da_method <- function(counts, method, group_1, group_2) {
  if (method == "propeller") return(run_sc_da_propeller_like(counts, group_1, group_2))
  if (method == "edger") return(run_sc_da_edger(counts, group_1, group_2))
  if (method == "chisq") return(run_sc_da_chisq(counts, group_1, group_2))
  run_sc_da_fisher(counts, group_1, group_2)
}

run_sc_da_fisher <- function(counts, group_1, group_2) {
  totals <- counts %>%
    group_by(.group) %>%
    summarise(total = sum(n_cells), .groups = "drop")
  
  ctypes <- sort(unique(counts$.celltype))
  
  res <- lapply(ctypes, function(ct) {
    g1_ct <- sum(counts$n_cells[counts$.group == group_1 & counts$.celltype == ct])
    g2_ct <- sum(counts$n_cells[counts$.group == group_2 & counts$.celltype == ct])
    
    g1_total <- totals$total[totals$.group == group_1]
    g2_total <- totals$total[totals$.group == group_2]
    
    mat <- matrix(
      c(g1_ct, g1_total - g1_ct, g2_ct, g2_total - g2_ct),
      nrow = 2,
      byrow = TRUE
    )
    
    ft <- fisher.test(mat)
    
    p1 <- (g1_ct + 0.5) / (g1_total + 1)
    p2 <- (g2_ct + 0.5) / (g2_total + 1)
    
    data.frame(
      cell_type = ct,
      group_1 = group_1,
      group_2 = group_2,
      group_1_count = g1_ct,
      group_2_count = g2_ct,
      group_1_prop = p1,
      group_2_prop = p2,
      logFC = log2(p2 / p1),
      P.Value = ft$p.value
    )
  })
  
  out <- bind_rows(res)
  out$FDR <- p.adjust(out$P.Value, method = "BH")
  out
}

run_sc_da_chisq <- function(counts, group_1, group_2) {
  totals <- counts %>%
    group_by(.group) %>%
    summarise(total = sum(n_cells), .groups = "drop")
  
  ctypes <- sort(unique(counts$.celltype))
  
  res <- lapply(ctypes, function(ct) {
    g1_ct <- sum(counts$n_cells[counts$.group == group_1 & counts$.celltype == ct])
    g2_ct <- sum(counts$n_cells[counts$.group == group_2 & counts$.celltype == ct])
    
    g1_total <- totals$total[totals$.group == group_1]
    g2_total <- totals$total[totals$.group == group_2]
    
    mat <- matrix(
      c(g1_ct, g1_total - g1_ct, g2_ct, g2_total - g2_ct),
      nrow = 2,
      byrow = TRUE
    )
    
    ctst <- suppressWarnings(chisq.test(mat))
    
    p1 <- (g1_ct + 0.5) / (g1_total + 1)
    p2 <- (g2_ct + 0.5) / (g2_total + 1)
    
    data.frame(
      cell_type = ct,
      group_1 = group_1,
      group_2 = group_2,
      group_1_count = g1_ct,
      group_2_count = g2_ct,
      group_1_prop = p1,
      group_2_prop = p2,
      logFC = log2(p2 / p1),
      P.Value = ctst$p.value
    )
  })
  
  out <- bind_rows(res)
  out$FDR <- p.adjust(out$P.Value, method = "BH")
  out
}

run_sc_da_propeller_like <- function(counts, group_1, group_2) {
  if (!requireNamespace("limma", quietly = TRUE)) {
    return(run_sc_da_fisher(counts, group_1, group_2))
  }
  
  prop_wide <- counts %>%
    select(.sample, .group, .celltype, proportion) %>%
    pivot_wider(
      names_from = .celltype,
      values_from = proportion,
      values_fill = 0
    ) %>%
    arrange(.sample)
  
  prop_mat <- as.matrix(prop_wide[, setdiff(colnames(prop_wide), c(".sample", ".group"))])
  rownames(prop_mat) <- prop_wide$.sample
  
  group <- factor(prop_wide$.group, levels = c(group_1, group_2))
  
  if (length(unique(group)) < 2) {
    return(run_sc_da_fisher(counts, group_1, group_2))
  }
  
  transformed <- asin(sqrt(prop_mat))
  design <- model.matrix(~ group)
  
  fit <- tryCatch(
    limma::lmFit(t(transformed), design),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    return(run_sc_da_fisher(counts, group_1, group_2))
  }
  
  fit <- limma::eBayes(fit)
  tab <- limma::topTable(fit, coef = 2, number = Inf, sort.by = "P")
  
  out <- data.frame(
    cell_type = rownames(tab),
    group_1 = group_1,
    group_2 = group_2,
    logFC = tab$logFC,
    P.Value = tab$P.Value,
    FDR = tab$adj.P.Val,
    row.names = NULL
  )
  
  add_group_summary(out, counts, group_1, group_2)
}

run_sc_da_edger <- function(counts, group_1, group_2) {
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    return(run_sc_da_fisher(counts, group_1, group_2))
  }
  
  count_wide <- counts %>%
    select(.sample, .group, .celltype, n_cells) %>%
    pivot_wider(
      names_from = .celltype,
      values_from = n_cells,
      values_fill = 0
    ) %>%
    arrange(.sample)
  
  sample_info <- count_wide %>% select(.sample, .group)
  
  count_mat <- as.matrix(count_wide[, setdiff(colnames(count_wide), c(".sample", ".group"))])
  rownames(count_mat) <- count_wide$.sample
  count_mat <- t(count_mat)
  
  group <- factor(sample_info$.group, levels = c(group_1, group_2))
  
  if (length(unique(group)) < 2) {
    return(run_sc_da_fisher(counts, group_1, group_2))
  }
  
  y <- edgeR::DGEList(counts = count_mat, group = group)
  y <- edgeR::calcNormFactors(y)
  
  design <- model.matrix(~ group)
  
  y <- tryCatch(edgeR::estimateDisp(y, design), error = function(e) NULL)
  if (is.null(y)) return(run_sc_da_fisher(counts, group_1, group_2))
  
  fit <- tryCatch(edgeR::glmQLFit(y, design), error = function(e) NULL)
  if (is.null(fit)) return(run_sc_da_fisher(counts, group_1, group_2))
  
  qlf <- edgeR::glmQLFTest(fit, coef = 2)
  tab <- edgeR::topTags(qlf, n = Inf)$table
  
  out <- data.frame(
    cell_type = rownames(tab),
    group_1 = group_1,
    group_2 = group_2,
    logFC = tab$logFC,
    P.Value = tab$PValue,
    FDR = tab$FDR,
    row.names = NULL
  )
  
  add_group_summary(out, counts, group_1, group_2)
}

add_group_summary <- function(out, counts, group_1, group_2) {
  agg <- counts %>%
    group_by(.group, .celltype) %>%
    summarise(n = sum(n_cells), .groups = "drop") %>%
    group_by(.group) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup()
  
  out %>%
    left_join(
      agg %>%
        filter(.group == group_1) %>%
        select(cell_type = .celltype, group_1_count = n, group_1_prop = prop),
      by = "cell_type"
    ) %>%
    left_join(
      agg %>%
        filter(.group == group_2) %>%
        select(cell_type = .celltype, group_2_count = n, group_2_prop = prop),
      by = "cell_type"
    )
}

# ==========================================================
# Plot helper functions
# ==========================================================

plot_sc_da_stacked <- function(long) {
  ggplot(
    long,
    aes(
      x = .sample,
      y = proportion,
      fill = .celltype,
      text = paste0(
        "Sample: ", .sample,
        "<br>Group: ", .group,
        "<br>Cell type: ", .celltype,
        "<br>Cells: ", n_cells,
        "<br>Proportion: ", round(proportion, 4)
      )
    )
  ) +
    geom_col(width = 0.85) +
    facet_grid(~ .group, scales = "free_x", space = "free_x") +
    scale_y_continuous(labels = percent_format()) +
    labs(
      x = "Sample",
      y = "Cell proportion",
      fill = "Cell type",
      title = "Cell type proportions per sample"
    ) +
    theme_classic(base_size = 13) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

plot_sc_da_boxplot <- function(long, selected_celltype) {
  dat <- long %>% filter(.celltype == selected_celltype)
  
  ggplot(
    dat,
    aes(
      x = .group,
      y = proportion,
      fill = .group,
      text = paste0(
        "Sample: ", .sample,
        "<br>Group: ", .group,
        "<br>Cell type: ", .celltype,
        "<br>Cells: ", n_cells,
        "<br>Proportion: ", round(proportion, 4)
      )
    )
  ) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.15, size = 2) +
    scale_y_continuous(labels = percent_format()) +
    labs(
      x = "Group",
      y = "Cell proportion",
      title = paste("Cell proportion:", selected_celltype)
    ) +
    theme_classic(base_size = 13) +
    theme(legend.position = "none")
}

plot_sc_da_volcano <- function(results, fdr_cutoff) {
  dat <- results %>%
    mutate(neg_log10_fdr = -log10(pmax(FDR, 1e-300)))
  
  ggplot(
    dat,
    aes(
      x = logFC,
      y = neg_log10_fdr,
      text = paste0(
        "Cell type: ", cell_type,
        "<br>logFC: ", round(logFC, 3),
        "<br>P value: ", signif(P.Value, 3),
        "<br>FDR: ", signif(FDR, 3),
        "<br>Significant: ", significant
      )
    )
  ) +
    geom_point(aes(shape = significant), size = 3, alpha = 0.85) +
    geom_hline(yintercept = -log10(fdr_cutoff), linetype = "dashed") +
    labs(
      x = "log2 fold change",
      y = "-log10 FDR",
      title = "Differential abundance volcano plot"
    ) +
    theme_classic(base_size = 13)
}

plot_sc_da_heatmap <- function(results) {
  dat <- results %>%
    mutate(cell_type = factor(cell_type, levels = cell_type[order(logFC)]))
  
  ggplot(
    dat,
    aes(
      x = "logFC",
      y = cell_type,
      fill = logFC,
      text = paste0(
        "Cell type: ", cell_type,
        "<br>logFC: ", round(logFC, 3),
        "<br>P value: ", signif(P.Value, 3),
        "<br>FDR: ", signif(FDR, 3)
      )
    )
  ) +
    geom_tile() +
    labs(
      x = "",
      y = "Cell type",
      fill = "logFC",
      title = "Differential abundance logFC heatmap"
    ) +
    theme_classic(base_size = 13) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
}

plot_sc_da_umap <- function(seurat_obj, metadata_col) {
  if (!"umap" %in% names(seurat_obj@reductions)) {
    stop("No UMAP reduction found.")
  }
  
  emb <- as.data.frame(Seurat::Embeddings(seurat_obj, reduction = "umap")[, 1:2, drop = FALSE])
  colnames(emb) <- c("UMAP_1", "UMAP_2")
  emb$metadata <- as.character(seurat_obj@meta.data[rownames(emb), metadata_col])
  
  ggplot(
    emb,
    aes(
      x = UMAP_1,
      y = UMAP_2,
      color = metadata,
      text = paste0(
        "Cell type: ", metadata,
        "<br>UMAP 1: ", round(UMAP_1, 3),
        "<br>UMAP 2: ", round(UMAP_2, 3)
      )
    )
  ) +
    geom_point(size = 0.8, alpha = 0.85) +
    labs(
      color = metadata_col,
      title = paste("UMAP colored by", metadata_col)
    ) +
    theme_classic(base_size = 13)
}
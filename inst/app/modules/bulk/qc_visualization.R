# ==========================================================
# modules/bulk/qc_visualization.R
# Bulk RNA-seq QC and visualization module
# - Interactive PCA
# - Sample correlation heatmap
# - Hierarchical clustering
# - Automatic plot updates when settings change
# - Functional PNG, PDF, and SVG downloads
# ==========================================================

mod_bulk_qc_ui <- function(id) {
  ns <- NS(id)

  tagList(
    shinyjs::useShinyjs(),

    actionButton(
      ns("help_toggle"),
      label = NULL,
      icon = icon("question-circle"),
      class = "btn btn-info btn-sm"
    ),

    shinyjs::hidden(
      div(
        id = ns("help_box"),
        style = paste(
          "border:1px solid #d9d9d9; padding:12px; margin-top:10px;",
          "border-radius:6px; background:#f7f7f7;"
        ),

        h4("Step-by-Step Guide"),

        h5("1. Choose a plot type"),
        p(paste(
          "Interactive PCA shows sample separation in reduced dimensions.",
          "The correlation heatmap shows similarity between samples.",
          "Hierarchical clustering shows sample relationships as a dendrogram."
        )),

        h5("2. Select the number of variable genes"),
        p(paste(
          "Choose how many highly variable genes to include.",
          "The plot updates automatically when this setting changes."
        )),

        h5("3. Review sample grouping"),
        p(paste(
          "PCA points and heatmap annotations use the experimental groups",
          "defined in the Groups module."
        )),

        h5("4. Explore the plot"),
        p(paste(
          "Hover over PCA points to view the sample name, group, and",
          "principal-component coordinates."
        )),

        h5("5. Download the plot"),
        p(paste(
          "Use the PNG, PDF, or SVG buttons to save the currently selected plot.",
          "PNG files are exported at 300 dpi."
        )),

        h5("6. Check sample and gene counts"),
        p(paste(
          "QC plots require at least two samples and at least two genes with",
          "non-zero variance after filtering."
        ))
      )
    ),

    br(),

    fluidRow(
      column(
        width = 4,
        selectInput(
          ns("plot_type"),
          "Select plot type",
          choices = c(
            "Interactive PCA",
            "Sample Correlation Heatmap",
            "Hierarchical Clustering"
          ),
          selected = "Interactive PCA"
        )
      ),
      column(
        width = 4,
        sliderInput(
          ns("n_top"),
          "Top variable genes",
          min = 500,
          max = 10000,
          value = 2000,
          step = 500
        )
      ),
      column(
        width = 4,
        tags$div(
          style = "padding-top:25px; display:flex; gap:8px; flex-wrap:wrap;",
          downloadButton(ns("dl_png"), "PNG", class = "btn btn-primary btn-sm"),
          downloadButton(ns("dl_pdf"), "PDF", class = "btn btn-primary btn-sm"),
          downloadButton(ns("dl_svg"), "SVG", class = "btn btn-primary btn-sm")
        )
      )
    ),

    hr(),
    uiOutput(ns("plot_ui")),
    br(),
    uiOutput(ns("interpretation_ui"))
  )
}

mod_bulk_qc_server <- function(id, bulk_data, groups) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    observeEvent(input$help_toggle, {
      shinyjs::toggle(id = "help_box", anim = TRUE)
    })

    observeEvent(input$interp_toggle, {
      shinyjs::toggle(id = "interpretation_body", anim = TRUE)
    })

    current_groups <- reactive({
      if (
        !is.null(groups) &&
        is.list(groups) &&
        !is.null(groups$groups) &&
        is.function(groups$groups)
      ) {
        out <- groups$groups()
        if (!is.null(out)) {
          return(out)
        }
      }

      data.frame(
        Group = character(0),
        Samples = character(0),
        stringsAsFactors = FALSE
      )
    })

    map_samples_to_groups <- function(groups_df, sample_names) {
      result <- setNames(rep("Unassigned", length(sample_names)), sample_names)

      if (
        is.null(groups_df) ||
        nrow(groups_df) == 0 ||
        !all(c("Group", "Samples") %in% colnames(groups_df))
      ) {
        return(result)
      }

      for (i in seq_len(nrow(groups_df))) {
        group_name <- trimws(as.character(groups_df$Group[i]))
        sample_text <- as.character(groups_df$Samples[i])

        if (
          is.na(group_name) ||
          !nzchar(group_name) ||
          is.na(sample_text) ||
          !nzchar(sample_text)
        ) {
          next
        }

        row_samples <- trimws(unlist(strsplit(sample_text, ";", fixed = TRUE)))
        row_samples <- row_samples[nzchar(row_samples)]
        matched <- intersect(row_samples, sample_names)

        if (length(matched) > 0) {
          result[matched] <- group_name
        }
      }

      result
    }

    analysis_matrix <- reactive({
      mat <- NULL

      if (
        !is.null(groups) &&
        is.list(groups) &&
        !is.null(groups$final) &&
        is.function(groups$final)
      ) {
        mat <- groups$final()
      }

      if (
        is.null(mat) &&
        !is.null(bulk_data) &&
        is.list(bulk_data) &&
        !is.null(bulk_data$final) &&
        is.function(bulk_data$final)
      ) {
        mat <- bulk_data$final()
      }

      req(mat)

      mat <- as.matrix(mat)
      suppressWarnings(storage.mode(mat) <- "numeric")

      validate(
        need(nrow(mat) >= 2, "QC visualization requires at least two genes."),
        need(ncol(mat) >= 2, paste(
          "QC visualization requires at least two samples after group filtering.",
          "Assign samples to groups or disable exclusion of ungrouped samples."
        )),
        need(!anyNA(mat), paste(
          "The expression matrix contains missing or non-numeric values.",
          "Please check the uploaded count matrix."
        )),
        need(all(is.finite(mat)), paste(
          "The expression matrix contains infinite values.",
          "Please check the uploaded count matrix."
        ))
      )

      mat
    })

    data_top <- reactive({
      mat <- analysis_matrix()
      log_mat <- log2(mat + 1)
      gene_variance <- apply(log_mat, 1, stats::var, na.rm = TRUE)

      valid <- which(is.finite(gene_variance) & gene_variance > 0)

      validate(
        need(length(valid) >= 2, paste(
          "At least two genes with non-zero variance are required.",
          "Check the uploaded data or sample filtering settings."
        ))
      )

      ordered <- valid[order(gene_variance[valid], decreasing = TRUE)]
      keep <- ordered[seq_len(min(input$n_top, length(ordered)))]

      log_mat[keep, , drop = FALSE]
    })

    pca_result <- reactive({
      m <- data_top()
      samples <- colnames(m)

      pca <- stats::prcomp(t(m), center = TRUE, scale. = TRUE)
      scores <- as.data.frame(pca$x, stringsAsFactors = FALSE)

      validate(
        need(all(c("PC1", "PC2") %in% colnames(scores)), paste(
          "PCA could not produce two principal components.",
          "At least two informative dimensions are required."
        ))
      )

      scores$Sample <- samples
      scores$Group <- factor(map_samples_to_groups(current_groups(), samples))

      variance <- pca$sdev^2
      variance_percent <- 100 * variance / sum(variance)

      list(
        pca = pca,
        scores = scores,
        variance_percent = variance_percent
      )
    })

    group_palette <- function(group_factor) {
      levels_present <- levels(droplevels(group_factor))
      n_groups <- length(levels_present)

      if (n_groups == 0) {
        return(character(0))
      }

      colours <- grDevices::hcl.colors(n_groups, palette = "Dark 3")
      stats::setNames(colours, levels_present)
    }

    pca_static_plot <- reactive({
      result <- pca_result()
      scores <- result$scores
      variance_percent <- result$variance_percent
      palette <- group_palette(scores$Group)

      ggplot2::ggplot(
        scores,
        ggplot2::aes(
          x = PC1,
          y = PC2,
          colour = Group,
          label = Sample
        )
      ) +
        ggplot2::geom_point(size = 3.6, alpha = 0.9) +
        ggplot2::geom_text(
          vjust = -0.9,
          size = 3.5,
          show.legend = FALSE,
          check_overlap = TRUE
        ) +
        ggplot2::scale_colour_manual(values = palette, drop = FALSE) +
        ggplot2::labs(
          title = "Principal Component Analysis",
          x = paste0("PC1 (", round(variance_percent[1], 1), "%)"),
          y = paste0("PC2 (", round(variance_percent[2], 1), "%)"),
          colour = "Group"
        ) +
        ggplot2::theme_minimal(base_size = 13) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
          legend.position = "right",
          panel.grid.minor = ggplot2::element_blank()
        )
    })

    heatmap_result <- reactive({
      m <- data_top()
      samples <- colnames(m)
      group_vector <- map_samples_to_groups(current_groups(), samples)

      annotation <- data.frame(
        Group = factor(group_vector),
        row.names = samples,
        check.names = FALSE
      )

      pheatmap::pheatmap(
        stats::cor(m, method = "pearson", use = "pairwise.complete.obs"),
        main = "Sample Correlation Heatmap",
        annotation_col = annotation,
        annotation_row = annotation,
        border_color = NA,
        silent = TRUE
      )
    })

    clustering_result <- reactive({
      stats::hclust(stats::dist(t(data_top()), method = "euclidean"))
    })

    output$plot_ui <- renderUI({
      analysis_matrix()

      if (identical(input$plot_type, "Interactive PCA")) {
        plotly::plotlyOutput(ns("pca_plot"), height = "550px")
      } else {
        plotOutput(ns("qc_plot_static"), height = "550px")
      }
    })

    output$pca_plot <- plotly::renderPlotly({
      result <- pca_result()
      scores <- result$scores
      variance_percent <- result$variance_percent

      tooltip_text <- paste0(
        "Sample: ", scores$Sample,
        "<br>Group: ", scores$Group,
        "<br>PC1: ", round(scores$PC1, 3),
        "<br>PC2: ", round(scores$PC2, 3)
      )

      palette <- group_palette(scores$Group)

      plotly::plot_ly(
        data = scores,
        x = ~PC1,
        y = ~PC2,
        color = ~Group,
        colors = unname(palette),
        text = tooltip_text,
        hoverinfo = "text",
        type = "scatter",
        mode = "markers+text",
        textposition = "top center",
        marker = list(
          size = 11,
          line = list(width = 1, color = "black")
        )
      ) |>
        plotly::layout(
          title = list(text = "Principal Component Analysis"),
          xaxis = list(
            title = paste0("PC1 (", round(variance_percent[1], 1), "%)")
          ),
          yaxis = list(
            title = paste0("PC2 (", round(variance_percent[2], 1), "%)")
          ),
          legend = list(title = list(text = "Group"))
        )
    })

    output$qc_plot_static <- renderPlot({
      if (identical(input$plot_type, "Sample Correlation Heatmap")) {
        heatmap <- heatmap_result()
        grid::grid.newpage()
        grid::grid.draw(heatmap$gtable)
      } else if (identical(input$plot_type, "Hierarchical Clustering")) {
        graphics::plot(
          clustering_result(),
          main = "Sample Hierarchical Clustering",
          xlab = "",
          sub = "",
          ylab = "Euclidean distance",
          hang = -1,
          cex = 0.9
        )
      }
    }, res = 110)

    interpretation_content <- reactive({
      req(input$plot_type)

      if (identical(input$plot_type, "Interactive PCA")) {
        tagList(
          h4("How to interpret PCA"),
          p("PCA shows whether samples separate according to the strongest expression patterns in the dataset."),
          h4("What to look for"),
          tags$ul(
            tags$li("Each point represents one sample."),
            tags$li("Samples close together have similar expression profiles."),
            tags$li("Samples far apart have more different expression profiles."),
            tags$li("Colours show the assigned experimental groups.")
          ),
          h4("Expected result"),
          p("Replicates from the same group should usually cluster together."),
          h4("Main caution"),
          p("A sample far from all other samples may indicate an outlier, batch effect, sample mix-up, or strong biological variation."),
          h4("Recommended next step"),
          p("Review unexpected sample separation before differential expression analysis.")
        )
      } else if (identical(input$plot_type, "Sample Correlation Heatmap")) {
        tagList(
          h4("How to interpret the correlation heatmap"),
          p("The heatmap shows similarity between samples based on their expression profiles."),
          h4("What to look for"),
          tags$ul(
            tags$li("Each cell shows the Pearson correlation between two samples."),
            tags$li("Higher correlation means greater similarity."),
            tags$li("Samples from the same group should usually show stronger correlation."),
            tags$li("Clear high-correlation blocks indicate consistent sample groups.")
          ),
          h4("Main caution"),
          p("A sample with low correlation to all other samples may be an outlier or poor-quality sample."),
          h4("Recommended next step"),
          p("Review unexpected low-correlation samples before downstream analysis.")
        )
      } else {
        tagList(
          h4("How to interpret hierarchical clustering"),
          p("The dendrogram shows sample relationships based on expression distance."),
          h4("What to look for"),
          tags$ul(
            tags$li("Samples joined by shorter branches are more similar."),
            tags$li("Samples joined by longer branches are more different."),
            tags$li("Replicates from the same group should ideally cluster together."),
            tags$li("Large separation can reflect biological or technical differences.")
          ),
          h4("Main caution"),
          p("Clustering by sequencing batch or another unwanted factor can indicate technical variation."),
          h4("Recommended next step"),
          p("Confirm that clustering agrees with the experimental design before downstream analysis.")
        )
      }
    })

    output$interpretation_ui <- renderUI({
      analysis_matrix()

      tags$div(
        style = paste(
          "border:1px solid #cfd8dc; border-radius:4px; margin-top:10px;",
          "background:white; overflow:hidden; box-shadow:0 1px 2px rgba(0,0,0,0.08);"
        ),
        tags$div(
          style = paste(
            "background:#1aa3b0; color:white; padding:9px 14px;",
            "display:flex; align-items:center; justify-content:space-between;"
          ),
          tags$span(
            icon("info-circle"),
            paste("Interpretation:", input$plot_type)
          ),
          actionButton(
            ns("interp_toggle"),
            label = NULL,
            icon = icon("minus"),
            class = "btn btn-xs",
            style = paste(
              "background:#147f8a; color:white; border:1px solid #0f6c75;",
              "padding:2px 8px;"
            )
          )
        ),
        div(
          id = ns("interpretation_body"),
          style = "padding:16px 18px; color:#111; background:white;",
          interpretation_content()
        )
      )
    })

    safe_plot_name <- reactive({
      switch(
        input$plot_type,
        "Interactive PCA" = "PCA",
        "Sample Correlation Heatmap" = "Correlation_Heatmap",
        "Hierarchical Clustering" = "Hierarchical_Clustering",
        "QC_Plot"
      )
    })

    draw_current_plot <- function() {
      if (identical(input$plot_type, "Interactive PCA")) {
        print(pca_static_plot())
      } else if (identical(input$plot_type, "Sample Correlation Heatmap")) {
        heatmap <- heatmap_result()
        grid::grid.newpage()
        grid::grid.draw(heatmap$gtable)
      } else if (identical(input$plot_type, "Hierarchical Clustering")) {
        graphics::plot(
          clustering_result(),
          main = "Sample Hierarchical Clustering",
          xlab = "",
          sub = "",
          ylab = "Euclidean distance",
          hang = -1,
          cex = 0.9
        )
      } else {
        stop("Unsupported QC plot type.")
      }
    }

    save_current_plot <- function(file, format) {
      data_top()

      format <- match.arg(format, c("png", "pdf", "svg"))

      if (identical(format, "png")) {
        grDevices::png(
          filename = file,
          width = 10,
          height = 8,
          units = "in",
          res = 300
        )
      } else if (identical(format, "pdf")) {
        grDevices::pdf(
          file = file,
          width = 10,
          height = 8,
          onefile = TRUE,
          useDingbats = FALSE
        )
      } else {
        if (requireNamespace("svglite", quietly = TRUE)) {
          svglite::svglite(file = file, width = 10, height = 8)
        } else {
          grDevices::svg(filename = file, width = 10, height = 8)
        }
      }

      on.exit(grDevices::dev.off(), add = TRUE)
      draw_current_plot()
    }

    output$dl_png <- downloadHandler(
      filename = function() {
        paste0(
          "CoTRA_Bulk_QC_",
          safe_plot_name(),
          "_",
          format(Sys.time(), "%Y%m%d_%H%M%S"),
          ".png"
        )
      },
      contentType = "image/png",
      content = function(file) {
        save_current_plot(file, "png")
      }
    )

    output$dl_pdf <- downloadHandler(
      filename = function() {
        paste0(
          "CoTRA_Bulk_QC_",
          safe_plot_name(),
          "_",
          format(Sys.time(), "%Y%m%d_%H%M%S"),
          ".pdf"
        )
      },
      contentType = "application/pdf",
      content = function(file) {
        save_current_plot(file, "pdf")
      }
    )

    output$dl_svg <- downloadHandler(
      filename = function() {
        paste0(
          "CoTRA_Bulk_QC_",
          safe_plot_name(),
          "_",
          format(Sys.time(), "%Y%m%d_%H%M%S"),
          ".svg"
        )
      },
      contentType = "image/svg+xml",
      content = function(file) {
        save_current_plot(file, "svg")
      }
    )
  })
}

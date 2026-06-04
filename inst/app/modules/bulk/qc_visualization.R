# ==========================================================
# modules/bulk/qc_visualization.R
# PCA, correlation heatmap, and dendrogram
# + interactive PCA and multi-format download
# + correct group-to-sample color mapping
# + help section
# + expandable plot interpretation section
# ==========================================================

mod_bulk_qc_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    shinyjs::useShinyjs(),
    
    actionButton(ns("help_toggle"), icon("question-circle"), class = "btn btn-info btn-sm"),
    
    shinyjs::hidden(
      div(
        id = ns("help_box"),
        style = "border:1px solid #d9d9d9; padding:12px; margin-top:10px; border-radius:6px; background:#f7f7f7;",
        
        h4("Step-by-Step Guide"),
        h5("1. Choose plot type"),
        p("Select one of the available QC plots. Interactive PCA shows sample clustering in reduced dimensions. Correlation heatmap shows similarity between samples. Hierarchical clustering shows sample relationships as a dendrogram."),
        
        h5("2. Select number of variable genes"),
        p("Choose how many top variable genes to use for the analysis. Higher values include more genes, while lower values focus on the most variable features."),
        
        h5("3. Redraw plot"),
        p("Click Redraw Plot to refresh the selected visualization using the current settings."),
        
        h5("4. Interpret sample grouping"),
        p("Samples are colored according to the assigned experimental groups. If ungrouped samples were excluded in the groups module, they will not appear here."),
        
        h5("5. Explore plot"),
        p("Interactive PCA allows hovering over samples to view sample name, group, and principal component coordinates. Static plots show heatmap or clustering structure."),
        
        h5("6. Download results"),
        p("Use the download button to save the current plot as PNG, PDF, or SVG for reports or publication figures."),
        
        h5("7. Check sample count"),
        p("QC plots require enough samples after filtering. If too few grouped samples remain, the module will show a message instead of a plot.")
      )
    ),
    
    br(),
    
    fluidRow(
      column(
        3,
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
        3,
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
        3,
        actionButton(ns("refresh"), "Redraw Plot", class = "btn btn-primary")
      ),
      column(
        3,
        dropdownButton(
          circle = TRUE,
          icon = icon("download"),
          tooltip = tooltipOptions(title = "Download"),
          downloadButton(ns("dl_png"), "PNG"),
          downloadButton(ns("dl_pdf"), "PDF"),
          downloadButton(ns("dl_svg"), "SVG")
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
    
    library(plotly)
    library(pheatmap)
    library(RColorBrewer)
    library(shinyjs)
    
    observeEvent(input$help_toggle, {
      shinyjs::toggle(id = "help_box", anim = TRUE)
    })
    
    observeEvent(input$interp_toggle, {
      shinyjs::toggle(id = "interpretation_body", anim = TRUE)
    })
    
    map_samples_to_groups <- function(groups_df, sample_names) {
      if (is.null(groups_df) || nrow(groups_df) == 0) {
        return(setNames(rep("Unassigned", length(sample_names)), sample_names))
      }
      
      all_samples <- unlist(strsplit(paste(groups_df$Samples, collapse = ";"), ";"))
      all_samples <- trimws(all_samples)
      
      counts_per_row <- sapply(strsplit(groups_df$Samples, ";"), length)
      all_groups <- rep(groups_df$Group, times = counts_per_row)
      
      mp <- setNames(all_groups, all_samples)
      
      out <- ifelse(sample_names %in% names(mp), mp[sample_names], "Unassigned")
      names(out) <- sample_names
      out
    }
    
    analysis_matrix <- reactive({
      if (!is.null(groups$final)) {
        mat <- groups$final()
      } else {
        mat <- bulk_data$final()
      }
      
      req(mat)
      mat
    })
    
    data_top <- reactive({
      mat <- analysis_matrix()
      
      req(ncol(mat) >= 2)
      req(nrow(mat) >= 2)
      
      m <- log2(mat + 1)
      vars <- apply(m, 1, var)
      top_idx <- order(vars, decreasing = TRUE)[1:min(input$n_top, nrow(m))]
      
      m[top_idx, , drop = FALSE]
    })
    
    output$plot_ui <- renderUI({
      mat <- analysis_matrix()
      
      if (ncol(mat) < 2) {
        return(
          tags$div(
            style = "padding:15px; color:#b22222;",
            "QC plots require at least 2 samples after group filtering. Please assign samples to groups or uncheck the exclude ungrouped option."
          )
        )
      }
      
      if (input$plot_type == "Interactive PCA") {
        plotly::plotlyOutput(ns("pca_plot"), height = "550px")
      } else {
        plotOutput(ns("qc_plot_static"), height = "550px")
      }
    })
    
    interpretation_content <- reactive({
      req(input$plot_type)
      
      if (input$plot_type == "Interactive PCA") {
        tagList(
          h4("How to interpret PCA"),
          p("PCA shows whether samples separate according to the strongest expression patterns in the dataset."),
          h4("What to look for"),
          tags$ul(
            tags$li("Each point represents one sample."),
            tags$li("Samples close together have similar expression profiles."),
            tags$li("Samples far apart are more different."),
            tags$li("Colors show the assigned experimental groups.")
          ),
          h4("Good result"),
          p("Replicates from the same group cluster together, and different groups show clear separation."),
          h4("Main caution"),
          p("A sample far away from all others may indicate an outlier, batch effect, sample mix-up, or strong biological variation."),
          h4("Recommended next step"),
          p("If groups separate clearly, continue to differential expression analysis. If unexpected patterns appear, check metadata, batch, and sample quality.")
        )
      } else if (input$plot_type == "Sample Correlation Heatmap") {
        tagList(
          h4("How to interpret the correlation heatmap"),
          p("The correlation heatmap shows how similar samples are based on their gene expression profiles."),
          h4("What to look for"),
          tags$ul(
            tags$li("Each cell shows the correlation between two samples."),
            tags$li("Higher correlation means samples are more similar."),
            tags$li("Samples from the same group should usually show stronger similarity."),
            tags$li("Clear blocks of high correlation suggest consistent biological groups.")
          ),
          h4("Good result"),
          p("Samples from the same group form strong correlation blocks."),
          h4("Main caution"),
          p("A sample with low correlation to all other samples may be an outlier or poor-quality sample."),
          h4("Recommended next step"),
          p("Review unexpected low-correlation samples before differential expression analysis.")
        )
      } else {
        tagList(
          h4("How to interpret hierarchical clustering"),
          p("The dendrogram shows sample relationships based on expression distance."),
          h4("What to look for"),
          tags$ul(
            tags$li("Samples joined by short branches are more similar."),
            tags$li("Samples joined by long branches are more different."),
            tags$li("Replicates from the same group should ideally cluster together."),
            tags$li("Large separation may reflect biological differences or technical effects.")
          ),
          h4("Good result"),
          p("Samples cluster according to the expected experimental design."),
          h4("Main caution"),
          p("If samples cluster by batch, sequencing run, or another unwanted factor, technical variation may affect the analysis."),
          h4("Recommended next step"),
          p("Check whether clustering agrees with the study design before downstream analysis.")
        )
      }
    })
    
    output$interpretation_ui <- renderUI({
      mat <- analysis_matrix()
      
      if (ncol(mat) < 2) {
        return(NULL)
      }
      
      tags$div(
        style = "border:1px solid #cfd8dc; border-radius:4px; margin-top:10px; background:white; overflow:hidden; box-shadow:0 1px 2px rgba(0,0,0,0.08);",
        
        tags$div(
          style = "background:#1aa3b0; color:white; padding:9px 14px; display:flex; align-items:center; justify-content:space-between;",
          
          tags$span(
            icon("info-circle"),
            paste("Interpretation:", input$plot_type)
          ),
          
          actionButton(
            ns("interp_toggle"),
            label = NULL,
            icon = icon("minus"),
            class = "btn btn-xs",
            style = "background:#147f8a; color:white; border:1px solid #0f6c75; padding:2px 8px;"
          )
        ),
        
        div(
          id = ns("interpretation_body"),
          style = "padding:16px 18px; color:#111; background:white;",
          interpretation_content()
        )
      )
    })
    
    output$pca_plot <- plotly::renderPlotly({
      req(data_top())
      
      m <- data_top()
      samples <- colnames(m)
      
      pca <- prcomp(t(m), scale. = TRUE)
      pcs <- as.data.frame(pca$x)
      pcs$Sample <- samples
      
      grp_vec <- map_samples_to_groups(groups$groups(), samples)
      pcs$Group <- factor(grp_vec)
      
      ve <- 100 * summary(pca)$importance[2, 1:2]
      pal <- brewer.pal(max(3, length(levels(pcs$Group))), "Set2")
      
      plotly::plot_ly(
        pcs,
        x = ~PC1,
        y = ~PC2,
        color = ~Group,
        colors = pal,
        text = ~paste(
          "Sample:", Sample,
          "<br>Group:", Group,
          "<br>PC1:", round(PC1, 2),
          "<br>PC2:", round(PC2, 2)
        ),
        hoverinfo = "text",
        type = "scatter",
        mode = "markers",
        marker = list(
          size = 12,
          line = list(width = 1, color = "black")
        )
      ) %>%
        plotly::layout(
          title = "Interactive PCA",
          xaxis = list(title = paste0("PC1 (", round(ve[1], 1), "%)")),
          yaxis = list(title = paste0("PC2 (", round(ve[2], 1), "%)"))
        )
    })
    
    output$qc_plot_static <- renderPlot({
      req(data_top())
      
      m <- data_top()
      samples <- colnames(m)
      grp_vec <- map_samples_to_groups(groups$groups(), samples)
      
      if (input$plot_type == "Sample Correlation Heatmap") {
        cor_mat <- cor(m)
        
        pheatmap(
          cor_mat,
          main = "Sample Correlation Heatmap",
          annotation_col = data.frame(
            Group = factor(grp_vec),
            row.names = samples
          )
        )
      } else if (input$plot_type == "Hierarchical Clustering") {
        d <- dist(t(m))
        hc <- hclust(d)
        
        plot(
          hc,
          main = "Sample Dendrogram",
          xlab = "",
          sub = ""
        )
      }
    })
    
    make_plot_for_download <- function() {
      m <- data_top()
      samples <- colnames(m)
      grp_vec <- map_samples_to_groups(groups$groups(), samples)
      
      if (input$plot_type == "Interactive PCA") {
        pca <- prcomp(t(m), scale. = TRUE)
        pcs <- as.data.frame(pca$x)
        pcs$Sample <- samples
        pcs$Group <- factor(grp_vec)
        
        pal <- brewer.pal(max(3, length(levels(pcs$Group))), "Set2")
        
        plot(
          pcs$PC1,
          pcs$PC2,
          pch = 19,
          col = pal[as.numeric(pcs$Group)],
          xlab = "PC1",
          ylab = "PC2",
          main = "PCA"
        )
        
        text(
          pcs$PC1,
          pcs$PC2,
          labels = pcs$Sample,
          pos = 3
        )
      } else if (input$plot_type == "Sample Correlation Heatmap") {
        cor_mat <- cor(m)
        
        pheatmap(
          cor_mat,
          main = "Sample Correlation Heatmap",
          annotation_col = data.frame(
            Group = factor(grp_vec),
            row.names = samples
          )
        )
      } else {
        d <- dist(t(m))
        hc <- hclust(d)
        
        plot(
          hc,
          main = "Sample Dendrogram",
          xlab = "",
          sub = ""
        )
      }
    }
    
    output$dl_png <- downloadHandler(
      filename = function() {
        "qc.png"
      },
      content = function(file) {
        png(file, width = 1200, height = 900, res = 150)
        make_plot_for_download()
        dev.off()
      }
    )
    
    output$dl_pdf <- downloadHandler(
      filename = function() {
        "qc.pdf"
      },
      content = function(file) {
        pdf(file, width = 10, height = 8)
        make_plot_for_download()
        dev.off()
      }
    )
    
    output$dl_svg <- downloadHandler(
      filename = function() {
        "qc.svg"
      },
      content = function(file) {
        svg(file, width = 10, height = 8)
        make_plot_for_download()
        dev.off()
      }
    )
  })
}
# ==========================================================
# modules/bulk/de_analysis.R
# Differential expression analysis using grouped samples only
# Supports DESeq2 and edgeR
# Includes help section, result downloads, PDF/SVG plot downloads
# Species selection added
# Comparison group auto updates when reference changes
# ==========================================================

suppressPackageStartupMessages({
  library(shiny)
  library(shinyjs)
  library(plotly)
  library(DT)
  library(htmlwidgets)
  library(webshot2)
})

mod_bulk_de_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    shinyjs::useShinyjs(),
    
    actionButton(ns("help_toggle"), icon("question-circle"), class = "btn btn-info btn-sm"),
    
    shinyjs::hidden(
      div(
        id = ns("help_box"),
        style = "border:1px solid #d9d9d9; padding:12px; margin-top:10px; border-radius:6px; background:#f7f7f7;",
        
        h4("Step-by-Step Guide"),
        
        h5("1. Define groups first"),
        p("Create experimental groups in the groups module before running differential expression analysis. Only samples assigned to groups are used here."),
        
        h5("2. Select analysis method"),
        p("Choose DESeq2 or edgeR. Both methods use raw count data and require grouped samples."),
        
        h5("2a. Select species"),
        p("Choose the organism for downstream compatibility with annotation and enrichment steps. Mouse is selected by default."),
        
        h5("3. Choose comparison groups"),
        p("Select one reference group and one comparison group. The analysis estimates differential expression for the comparison group relative to the reference group."),
        
        h5("4. Run analysis"),
        p("Click Run DE analysis to calculate log2 fold changes, adjusted p values, and significance statistics for all genes."),
        
        h5("5. Adjust significance thresholds"),
        p("Use the adjusted p value cutoff and absolute log2 fold change cutoff to filter significant genes."),
        
        h5("6. Explore results"),
        p("The volcano plot, MA plot, PCA, heatmap, and result tables help inspect significant genes and sample separation between selected groups."),
        
        h5("7. Search and download"),
        p("Use the gene search box to filter by gene name. Download tables or export plots as PDF or SVG.")
      )
    ),
    
    br(),
    
    fluidRow(
      column(
        3,
        selectInput(
          ns("species"),
          "Species",
          choices = c("Mouse" = "mouse", "Human" = "human", "Rat" = "rat"),
          selected = "mouse"
        )
      ),
      column(
        3,
        selectInput(
          ns("method"),
          "DE method",
          choices = c("DESeq2", "edgeR"),
          selected = "DESeq2"
        )
      ),
      column(3, uiOutput(ns("ref_group_ui"))),
      column(3, uiOutput(ns("comp_group_ui")))
    ),
    
    br(),
    
    fluidRow(
      column(
        3,
        actionButton(
          ns("run_de"),
          "Run DE analysis",
          class = "btn btn-primary"
        )
      ),
      column(
        3,
        sliderInput(
          ns("p_cut_slider"),
          "Adjusted p cutoff",
          min = 0,
          max = 1,
          value = 0.05,
          step = 0.005
        ),
        numericInput(
          ns("p_cut_num"),
          NULL,
          value = 0.05,
          min = 0,
          max = 1,
          step = 0.001
        )
      ),
      column(
        3,
        sliderInput(
          ns("fc_cut_slider"),
          "Absolute log2FC cutoff",
          min = 0,
          max = 5,
          value = 1,
          step = 0.1
        ),
        numericInput(
          ns("fc_cut_num"),
          NULL,
          value = 1,
          min = 0,
          max = 5,
          step = 0.01
        )
      ),
      column(
        3,
        textInput(
          ns("gene_search"),
          "Search gene",
          value = ""
        ),
        br(),
        verbatimTextOutput(ns("status"))
      )
    ),
    
    hr(),
    
    fluidRow(
      column(
        4,
        downloadButton(
          ns("download_csv"),
          "Download filtered DE (CSV)",
          class = "btn btn-success"
        )
      ),
      column(
        4,
        downloadButton(
          ns("download_all_csv"),
          "Download all DE results (CSV)",
          class = "btn btn-default"
        )
      )
    ),
    
    hr(),
    
    tabsetPanel(
      tabPanel(
        "Volcano Plot",
        br(),
        fluidRow(
          column(3, downloadButton(ns("download_volcano_pdf"), "Download PDF")),
          column(3, downloadButton(ns("download_volcano_svg"), "Download SVG"))
        ),
        br(),
        plotly::plotlyOutput(ns("volcano"), height = "550px")
      ),
      
      tabPanel(
        "MA Plot",
        br(),
        fluidRow(
          column(3, downloadButton(ns("download_ma_pdf"), "Download PDF")),
          column(3, downloadButton(ns("download_ma_svg"), "Download SVG"))
        ),
        br(),
        plotly::plotlyOutput(ns("ma_plot"), height = "550px")
      ),
      
      tabPanel(
        "PCA of significant genes",
        br(),
        fluidRow(
          column(3, downloadButton(ns("download_pca_pdf"), "Download PDF")),
          column(3, downloadButton(ns("download_pca_svg"), "Download SVG"))
        ),
        br(),
        plotly::plotlyOutput(ns("pca_de"), height = "550px")
      ),
      
      tabPanel(
        "Heatmap of significant genes",
        br(),
        numericInput(
          ns("n_genes_heat"),
          "Max genes in heatmap",
          value = 50,
          min = 10,
          max = 200,
          step = 10
        ),
        fluidRow(
          column(3, downloadButton(ns("download_heatmap_pdf"), "Download PDF")),
          column(3, downloadButton(ns("download_heatmap_svg"), "Download SVG"))
        ),
        br(),
        plotly::plotlyOutput(ns("heatmap"), height = "700px")
      ),
      
      tabPanel(
        "DE Tables",
        tabsetPanel(
          tabPanel("All filtered", DT::DTOutput(ns("table_sig"))),
          tabPanel("Upregulated", DT::DTOutput(ns("table_up"))),
          tabPanel("Downregulated", DT::DTOutput(ns("table_down"))),
          tabPanel("All genes", DT::DTOutput(ns("table_all")))
        )
      )
    )
  )
}

mod_bulk_de_server <- function(id, bulk_data, groups) {
  moduleServer(id, function(input, output, session) {
    
    ns <- session$ns
    
    rv <- reactiveValues(
      results = NULL,
      metadata = NULL
    )
    
    observeEvent(input$help_toggle, {
      shinyjs::toggle(id = "help_box", anim = TRUE)
    })
    
    apply_gene_search <- function(df, query) {
      if (is.null(df)) return(NULL)
      if (is.null(query) || trimws(query) == "") return(df)
      
      idx <- grepl(trimws(query), df$gene, ignore.case = TRUE)
      df[idx, , drop = FALSE]
    }
    
    safe_filename <- function(x) {
      x <- gsub("[^A-Za-z0-9_\\-]+", "_", x)
      x <- gsub("_+", "_", x)
      x
    }
    
    export_plotly_widget <- function(p, file, width = 1100, height = 800) {
      tmp_html <- tempfile(fileext = ".html")
      
      htmlwidgets::saveWidget(
        widget = p,
        file = tmp_html,
        selfcontained = TRUE
      )
      
      webshot2::webshot(
        url = tmp_html,
        file = file,
        vwidth = width,
        vheight = height
      )
    }
    
    p_cut <- reactive({
      v <- input$p_cut_num
      if (is.null(v) || is.na(v)) return(0.05)
      max(0, min(1, v))
    })
    
    fc_cut <- reactive({
      v <- input$fc_cut_num
      if (is.null(v) || is.na(v)) return(1)
      max(0, min(5, v))
    })
    
    observeEvent(input$p_cut_slider, {
      updateNumericInput(session, "p_cut_num", value = input$p_cut_slider)
    }, ignoreInit = TRUE)
    
    observeEvent(input$p_cut_num, {
      updateSliderInput(session, "p_cut_slider", value = p_cut())
    }, ignoreInit = TRUE)
    
    observeEvent(input$fc_cut_slider, {
      updateNumericInput(session, "fc_cut_num", value = input$fc_cut_slider)
    }, ignoreInit = TRUE)
    
    observeEvent(input$fc_cut_num, {
      updateSliderInput(session, "fc_cut_slider", value = fc_cut())
    }, ignoreInit = TRUE)
    
    groups_df <- reactive({
      req(groups$groups())
      groups$groups()
    })
    
    available_groups <- reactive({
      df <- groups_df()
      
      if (is.null(df) || nrow(df) == 0) {
        return(character(0))
      }
      
      unique(trimws(df$Group))
    })
    
    output$ref_group_ui <- renderUI({
      grps <- available_groups()
      
      selectInput(
        ns("ref_group"),
        "Reference group",
        choices = grps,
        selected = if (length(grps) >= 1) grps[1] else NULL
      )
    })
    
    output$comp_group_ui <- renderUI({
      grps <- available_groups()
      default_comp <- NULL
      
      if (length(grps) >= 2) {
        if (!is.null(input$ref_group) && input$ref_group %in% grps) {
          other_grps <- setdiff(grps, input$ref_group)
          default_comp <- other_grps[1]
        } else {
          default_comp <- grps[2]
        }
      }
      
      selectInput(
        ns("comp_group"),
        "Comparison group",
        choices = grps,
        selected = default_comp
      )
    })
    
    observeEvent(input$ref_group, {
      grps <- available_groups()
      req(length(grps) >= 2)
      
      other_grps <- setdiff(grps, input$ref_group)
      new_comp <- if (length(other_grps) > 0) other_grps[1] else NULL
      
      updateSelectInput(
        session,
        "comp_group",
        choices = grps,
        selected = new_comp
      )
    }, ignoreInit = TRUE)
    
    grouped_sample_map <- reactive({
      df <- groups_df()
      req(nrow(df) > 0)
      
      out <- do.call(
        rbind,
        lapply(seq_len(nrow(df)), function(i) {
          smp <- trimws(unlist(strsplit(df$Samples[i], ";")))
          smp <- smp[nzchar(smp)]
          
          if (length(smp) == 0) {
            return(NULL)
          }
          
          data.frame(
            Sample = smp,
            Group = trimws(df$Group[i]),
            stringsAsFactors = FALSE
          )
        })
      )
      
      out <- out[!duplicated(out$Sample), , drop = FALSE]
      rownames(out) <- out$Sample
      out
    })
    
    grouped_norm_matrix <- reactive({
      mat <- bulk_data$final()
      req(mat)
      
      smap <- grouped_sample_map()
      keep <- intersect(colnames(mat), smap$Sample)
      
      validate(
        need(length(keep) >= 2, "At least two grouped samples are required for DE analysis.")
      )
      
      mat[, keep, drop = FALSE]
    })
    
    grouped_raw_matrix <- reactive({
      mat <- bulk_data$raw()
      req(mat)
      
      smap <- grouped_sample_map()
      keep <- intersect(colnames(mat), smap$Sample)
      
      validate(
        need(length(keep) >= 2, "At least two grouped samples are required for DE analysis.")
      )
      
      mat[, keep, drop = FALSE]
    })
    
    comparison_data <- reactive({
      req(input$ref_group, input$comp_group)
      
      validate(
        need(input$ref_group != input$comp_group, "Choose two different groups for comparison.")
      )
      
      smap <- grouped_sample_map()
      
      sel <- smap$Group %in% c(input$ref_group, input$comp_group)
      smap2 <- smap[sel, , drop = FALSE]
      
      validate(
        need(nrow(smap2) >= 2, "Selected groups do not contain enough samples."),
        need(sum(smap2$Group == input$ref_group) >= 1, "Reference group has no samples."),
        need(sum(smap2$Group == input$comp_group) >= 1, "Comparison group has no samples.")
      )
      
      raw_mat <- grouped_raw_matrix()[, smap2$Sample, drop = FALSE]
      norm_mat <- grouped_norm_matrix()[, smap2$Sample, drop = FALSE]
      
      smap2$Group <- factor(
        smap2$Group,
        levels = c(input$ref_group, input$comp_group)
      )
      
      list(
        raw = raw_mat,
        norm = norm_mat,
        sample_info = smap2
      )
    })
    
    observeEvent(input$run_de, {
      req(comparison_data())
      
      dat <- comparison_data()
      
      raw_mat <- dat$raw
      norm_mat <- dat$norm
      sample_info <- dat$sample_info
      
      validate(
        need(ncol(raw_mat) >= 2, "At least two grouped samples are required."),
        need(length(unique(sample_info$Group)) == 2, "Exactly two groups are required for DE analysis."),
        need(all(sample_info$Sample == colnames(raw_mat)), "Sample order mismatch detected."),
        need(sum(sample_info$Group == input$ref_group) >= 1, "Reference group has no samples."),
        need(sum(sample_info$Group == input$comp_group) >= 1, "Comparison group has no samples.")
      )
      
      output$status <- renderText("Running differential expression analysis...")
      
      res_df <- NULL
      
      if (identical(input$method, "DESeq2")) {
        coldata <- data.frame(group = sample_info$Group)
        rownames(coldata) <- sample_info$Sample
        
        dds <- DESeq2::DESeqDataSetFromMatrix(
          countData = round(as.matrix(raw_mat)),
          colData = coldata,
          design = ~ group
        )
        
        keep <- rowSums(DESeq2::counts(dds) >= 10) >= 2
        
        validate(
          need(sum(keep) >= 2, "Too few genes passed count filtering for DESeq2.")
        )
        
        dds <- dds[keep, ]
        dds <- DESeq2::DESeq(dds)
        
        res <- DESeq2::results(
          dds,
          contrast = c("group", input$comp_group, input$ref_group)
        )
        
        res_df <- as.data.frame(res)
        res_df$gene <- rownames(res_df)
        res_df$baseMean <- res_df$baseMean
        res_df$lfc <- res_df$log2FoldChange
        res_df$pvalue <- res_df$pvalue
        res_df$padj <- res_df$padj
      }
      
      if (identical(input$method, "edgeR")) {
        y <- edgeR::DGEList(
          counts = round(as.matrix(raw_mat)),
          group = sample_info$Group
        )
        
        keep <- edgeR::filterByExpr(y, group = sample_info$Group)
        
        validate(
          need(sum(keep) >= 2, "Too few genes passed count filtering for edgeR.")
        )
        
        y <- y[keep, , keep.lib.sizes = FALSE]
        y <- edgeR::calcNormFactors(y)
        
        design <- model.matrix(
          ~ group,
          data = data.frame(group = sample_info$Group)
        )
        
        y <- edgeR::estimateDisp(y, design)
        fit <- edgeR::glmQLFit(y, design)
        test <- edgeR::glmQLFTest(fit, coef = 2)
        tab <- edgeR::topTags(test, n = Inf)$table
        
        res_df <- as.data.frame(tab)
        res_df$gene <- rownames(res_df)
        res_df$baseMean <- rowMeans(edgeR::cpm(y, log = FALSE))
        res_df$lfc <- res_df$logFC
        res_df$pvalue <- res_df$PValue
        res_df$padj <- res_df$FDR
      }
      
      validate(
        need(!is.null(res_df), "DE analysis failed to produce results.")
      )
      
      res_df$padj[is.na(res_df$padj)] <- 1
      res_df$pvalue[is.na(res_df$pvalue)] <- 1
      res_df$lfc[is.na(res_df$lfc)] <- 0
      
      res_df$significant <- abs(res_df$lfc) >= fc_cut() & res_df$padj <= p_cut()
      
      res_df <- res_df[order(res_df$padj, -abs(res_df$lfc)), , drop = FALSE]
      
      rv$results <- res_df
      
      rv$metadata <- list(
        species = input$species,
        method = input$method,
        reference_group = input$ref_group,
        comparison_group = input$comp_group,
        sample_info = sample_info,
        norm_matrix = norm_mat
      )
      
      output$status <- renderText({
        paste0(
          "DE complete. ",
          nrow(res_df), " genes tested. ",
          sum(abs(res_df$lfc) >= fc_cut() & res_df$padj <= p_cut(), na.rm = TRUE),
          " genes pass current thresholds."
        )
      })
    })
    
    de_full <- reactive({
      rv$results
    })
    
    de_sig <- reactive({
      df <- de_full()
      
      if (is.null(df)) {
        return(NULL)
      }
      
      df <- df[abs(df$lfc) >= fc_cut() & df$padj <= p_cut(), , drop = FALSE]
      df
    })
    
    de_sig_search <- reactive({
      apply_gene_search(de_sig(), input$gene_search)
    })
    
    de_full_search <- reactive({
      apply_gene_search(de_full(), input$gene_search)
    })
    
    de_up <- reactive({
      df <- de_sig_search()
      
      if (is.null(df)) {
        return(NULL)
      }
      
      df[df$lfc > 0, , drop = FALSE]
    })
    
    de_down <- reactive({
      df <- de_sig_search()
      
      if (is.null(df)) {
        return(NULL)
      }
      
      df[df$lfc < 0, , drop = FALSE]
    })
    
    output$table_sig <- DT::renderDT({
      DT::datatable(
        de_sig_search(),
        rownames = FALSE,
        options = list(scrollX = TRUE, pageLength = 15)
      )
    })
    
    output$table_up <- DT::renderDT({
      DT::datatable(
        de_up(),
        rownames = FALSE,
        options = list(scrollX = TRUE, pageLength = 15)
      )
    })
    
    output$table_down <- DT::renderDT({
      DT::datatable(
        de_down(),
        rownames = FALSE,
        options = list(scrollX = TRUE, pageLength = 15)
      )
    })
    
    output$table_all <- DT::renderDT({
      DT::datatable(
        de_full_search(),
        rownames = FALSE,
        options = list(scrollX = TRUE, pageLength = 15)
      )
    })
    
    volcano_plot_obj <- reactive({
      req(de_full())
      
      df <- de_full_search()
      
      validate(
        need(!is.null(df), "No DE results available."),
        need(nrow(df) > 0, "No genes available for plotting.")
      )
      
      df$plot_group <- ifelse(
        df$padj <= p_cut() & abs(df$lfc) >= fc_cut(),
        "Significant",
        "Not significant"
      )
      
      df$neglog10padj <- -log10(df$padj + 1e-300)
      
      plotly::plot_ly(
        df,
        x = ~lfc,
        y = ~neglog10padj,
        type = "scatter",
        mode = "markers",
        color = ~plot_group,
        colors = c("Not significant" = "grey", "Significant" = "red"),
        text = ~paste0(
          "Gene: ", gene,
          "<br>log2FC: ", round(lfc, 3),
          "<br>p value: ", signif(pvalue, 3),
          "<br>Adj p: ", signif(padj, 3)
        ),
        hoverinfo = "text"
      ) %>%
        plotly::layout(
          title = paste0("Volcano Plot: ", input$comp_group, " vs ", input$ref_group),
          xaxis = list(title = "log2 fold change"),
          yaxis = list(title = "-log10 adjusted p value")
        )
    })
    
    ma_plot_obj <- reactive({
      req(de_full())
      
      df <- de_full_search()
      
      validate(
        need(!is.null(df), "No DE results available."),
        need(nrow(df) > 0, "No genes available for plotting.")
      )
      
      df$plot_group <- ifelse(
        df$padj <= p_cut() & abs(df$lfc) >= fc_cut(),
        "Significant",
        "Not significant"
      )
      
      plotly::plot_ly(
        df,
        x = ~log10(baseMean + 1),
        y = ~lfc,
        type = "scatter",
        mode = "markers",
        color = ~plot_group,
        colors = c("Not significant" = "grey", "Significant" = "red"),
        text = ~paste0(
          "Gene: ", gene,
          "<br>Mean expression: ", round(baseMean, 3),
          "<br>log2FC: ", round(lfc, 3),
          "<br>p value: ", signif(pvalue, 3),
          "<br>Adj p: ", signif(padj, 3)
        ),
        hoverinfo = "text"
      ) %>%
        plotly::layout(
          title = paste0("MA Plot: ", input$comp_group, " vs ", input$ref_group),
          xaxis = list(title = "log10 mean expression"),
          yaxis = list(title = "log2 fold change")
        )
    })
    
    pca_plot_obj <- reactive({
      req(rv$metadata, de_sig())
      
      sig <- de_sig()
      md <- rv$metadata
      
      validate(
        need(!is.null(sig), "No significant genes available."),
        need(nrow(sig) >= 2, "At least two significant genes are needed for PCA.")
      )
      
      mat <- md$norm_matrix
      genes <- intersect(sig$gene, rownames(mat))
      
      validate(
        need(length(genes) >= 2, "Significant genes are not available in the expression matrix.")
      )
      
      genes <- genes[seq_len(min(length(genes), 500))]
      
      sub <- mat[genes, , drop = FALSE]
      sub <- log2(sub + 1)
      sub <- t(scale(t(sub)))
      sub[is.na(sub)] <- 0
      sub[is.infinite(sub)] <- 0
      
      pca <- prcomp(t(sub), scale. = FALSE)
      
      pcs <- as.data.frame(pca$x)
      pcs$Sample <- rownames(pcs)
      pcs$Group <- md$sample_info[pcs$Sample, "Group"]
      
      ve <- 100 * summary(pca)$importance[2, 1:2]
      
      plotly::plot_ly(
        pcs,
        x = ~PC1,
        y = ~PC2,
        type = "scatter",
        mode = "markers+text",
        color = ~Group,
        text = ~Sample,
        textposition = "top center",
        hoverinfo = "text",
        hovertext = ~paste0(
          "Sample: ", Sample,
          "<br>Group: ", Group
        )
      ) %>%
        plotly::layout(
          title = "PCA of significant genes",
          xaxis = list(title = paste0("PC1 (", round(ve[1], 1), "%)")),
          yaxis = list(title = paste0("PC2 (", round(ve[2], 1), "%)"))
        )
    })
    
    heatmap_plot_obj <- reactive({
      req(rv$metadata, de_sig())
      
      sig <- de_sig()
      md <- rv$metadata
      
      validate(
        need(!is.null(sig), "No significant genes available."),
        need(nrow(sig) > 0, "No significant genes available for heatmap.")
      )
      
      mat <- md$norm_matrix
      genes <- intersect(sig$gene, rownames(mat))
      
      validate(
        need(length(genes) > 0, "Significant genes are not available in the expression matrix.")
      )
      
      n_heat <- input$n_genes_heat
      
      if (is.null(n_heat) || is.na(n_heat)) {
        n_heat <- 50
      }
      
      genes <- genes[seq_len(min(length(genes), n_heat))]
      
      sub <- mat[genes, , drop = FALSE]
      sub <- log2(sub + 1)
      sub <- t(scale(t(sub)))
      sub[is.na(sub)] <- 0
      sub[is.infinite(sub)] <- 0
      
      if (nrow(sub) > 1) {
        row_ord <- hclust(dist(sub))$order
        sub <- sub[row_ord, , drop = FALSE]
      }
      
      if (ncol(sub) > 1) {
        col_ord <- hclust(dist(t(sub)))$order
        sub <- sub[, col_ord, drop = FALSE]
      }
      
      plotly::plot_ly(
        x = colnames(sub),
        y = rownames(sub),
        z = sub,
        type = "heatmap",
        colorscale = "RdBu",
        reversescale = TRUE
      ) %>%
        plotly::layout(
          title = paste0("Heatmap of significant genes (", nrow(sub), " genes)"),
          xaxis = list(title = "Samples"),
          yaxis = list(title = "Genes", autorange = "reversed")
        )
    })
    
    output$volcano <- plotly::renderPlotly({
      volcano_plot_obj()
    })
    
    output$ma_plot <- plotly::renderPlotly({
      ma_plot_obj()
    })
    
    output$pca_de <- plotly::renderPlotly({
      pca_plot_obj()
    })
    
    output$heatmap <- plotly::renderPlotly({
      heatmap_plot_obj()
    })
    
    output$download_csv <- downloadHandler(
      filename = function() {
        paste0(
          "DE_filtered_",
          safe_filename(input$comp_group),
          "_vs_",
          safe_filename(input$ref_group),
          "_",
          safe_filename(input$method),
          ".csv"
        )
      },
      content = function(file) {
        write.csv(de_sig_search(), file, row.names = FALSE)
      }
    )
    
    output$download_all_csv <- downloadHandler(
      filename = function() {
        paste0(
          "DE_all_",
          safe_filename(input$comp_group),
          "_vs_",
          safe_filename(input$ref_group),
          "_",
          safe_filename(input$method),
          ".csv"
        )
      },
      content = function(file) {
        write.csv(de_full_search(), file, row.names = FALSE)
      }
    )
    
    output$download_volcano_pdf <- downloadHandler(
      filename = function() {
        paste0(
          "Volcano_",
          safe_filename(input$comp_group),
          "_vs_",
          safe_filename(input$ref_group),
          ".pdf"
        )
      },
      content = function(file) {
        export_plotly_widget(
          volcano_plot_obj(),
          file,
          width = 1100,
          height = 800
        )
      }
    )
    
    output$download_volcano_svg <- downloadHandler(
      filename = function() {
        paste0(
          "Volcano_",
          safe_filename(input$comp_group),
          "_vs_",
          safe_filename(input$ref_group),
          ".svg"
        )
      },
      content = function(file) {
        export_plotly_widget(
          volcano_plot_obj(),
          file,
          width = 1100,
          height = 800
        )
      }
    )
    
    output$download_ma_pdf <- downloadHandler(
      filename = function() {
        paste0(
          "MA_plot_",
          safe_filename(input$comp_group),
          "_vs_",
          safe_filename(input$ref_group),
          ".pdf"
        )
      },
      content = function(file) {
        export_plotly_widget(
          ma_plot_obj(),
          file,
          width = 1100,
          height = 800
        )
      }
    )
    
    output$download_ma_svg <- downloadHandler(
      filename = function() {
        paste0(
          "MA_plot_",
          safe_filename(input$comp_group),
          "_vs_",
          safe_filename(input$ref_group),
          ".svg"
        )
      },
      content = function(file) {
        export_plotly_widget(
          ma_plot_obj(),
          file,
          width = 1100,
          height = 800
        )
      }
    )
    
    output$download_pca_pdf <- downloadHandler(
      filename = function() {
        paste0(
          "PCA_significant_genes_",
          safe_filename(input$comp_group),
          "_vs_",
          safe_filename(input$ref_group),
          ".pdf"
        )
      },
      content = function(file) {
        export_plotly_widget(
          pca_plot_obj(),
          file,
          width = 1100,
          height = 800
        )
      }
    )
    
    output$download_pca_svg <- downloadHandler(
      filename = function() {
        paste0(
          "PCA_significant_genes_",
          safe_filename(input$comp_group),
          "_vs_",
          safe_filename(input$ref_group),
          ".svg"
        )
      },
      content = function(file) {
        export_plotly_widget(
          pca_plot_obj(),
          file,
          width = 1100,
          height = 800
        )
      }
    )
    
    output$download_heatmap_pdf <- downloadHandler(
      filename = function() {
        paste0(
          "Heatmap_significant_genes_",
          safe_filename(input$comp_group),
          "_vs_",
          safe_filename(input$ref_group),
          ".pdf"
        )
      },
      content = function(file) {
        export_plotly_widget(
          heatmap_plot_obj(),
          file,
          width = 1200,
          height = 1000
        )
      }
    )
    
    output$download_heatmap_svg <- downloadHandler(
      filename = function() {
        paste0(
          "Heatmap_significant_genes_",
          safe_filename(input$comp_group),
          "_vs_",
          safe_filename(input$ref_group),
          ".svg"
        )
      },
      content = function(file) {
        export_plotly_widget(
          heatmap_plot_obj(),
          file,
          width = 1200,
          height = 1000
        )
      }
    )
    
    return(reactive(rv$results))
  })
}
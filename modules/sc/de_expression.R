mod_sc_de_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    fluidRow(
      h2("To start the differential expression analysis, first select the data type and then define the groups of your interest"),
      hr(),
      
      box(
        title = p("Parameters to Define the Groups of Cells",
                  actionButton(ns("HELP_DEG_INPUT"), "", icon = icon("question-circle"), class = "btn-xs")),
        width = 3,
        height = NULL,
        solidHeader = TRUE,
        status = "primary",
        
        uiOutput(ns("select_meta_fac1")),
        uiOutput(ns("range_meta_fac1")),
        
        uiOutput(ns("select_meta_fac2")),
        uiOutput(ns("range_meta_fac2")),
        
        uiOutput(ns("select_reduc")),
        uiOutput(ns("dea_vis_btn")),
        hr()
      ),
      
      conditionalPanel(
        condition = "input.dea_vis_btn2",
        ns = ns,
        
        box(
          title = p("Parameters for Differential Expression Analysis",
                    actionButton(ns("HELP_DEG_PARAMS"), "", icon = icon("question-circle"), class = "btn-xs")),
          width = 3,
          height = NULL,
          solidHeader = TRUE,
          status = "primary",
          
          uiOutput(ns("min_pct")),
          uiOutput(ns("nb_logfc")),
          uiOutput(ns("select_method_Markers")),
          uiOutput(ns("checkbox_only_pos")),
          useShinyalert(),
          uiOutput(ns("dea_run_btn"))
        ),
        
        box(
          title = p("Visualization of the Groups of Cells",
                    actionButton(ns("HELP_DEG_VIS"), "", icon = icon("question-circle"), class = "btn-xs")),
          width = 6,
          height = NULL,
          solidHeader = TRUE,
          status = "primary",
          
          plotlyOutput(ns("dea_vis_plot"), height = 450),
          br(),
          div(
            style = "padding-left: 8em; font-size: 12pt; color: black;",
            textOutput(ns("state_vis_out"))
          ),
          br(),
          radioButtons(
            ns("fmt_dea_vis"),
            "Download format",
            choices = c("PDF" = "pdf", "SVG" = "svg"),
            inline = TRUE
          ),
          downloadButton(ns("down_dea_vis"), "Download group plot", class = "btn btn-primary btn-sm")
        )
      )
    ),
    
    fluidRow(
      br(),
      br(),
      
      conditionalPanel(
        condition = "input.dea_vis_btn2",
        ns = ns,
        
        conditionalPanel(
          condition = "input.dea_run_btn2",
          ns = ns,
          
          box(
            title = p("Differentially Expressed Genes",
                      actionButton(ns("HELP_DEG_TABLE"), "", icon = icon("question-circle"), class = "btn-xs")),
            height = NULL,
            solidHeader = TRUE,
            status = "primary",
            width = 14,
            
            tabsetPanel(
              id = ns("panel_deg"),
              
              tabPanel(
                "Table of Differentially Expressed Genes",
                box(
                  width = NULL,
                  height = NULL,
                  solidHeader = TRUE,
                  status = "primary",
                  br(),
                  div(
                    style = "padding-left: 6em; font-size: 16pt; color: black;",
                    htmlOutput(ns("dea_table_stat"))
                  ),
                  br(),
                  div(
                    style = "padding-left: 6em; font-size: 12pt; color: black;",
                    textOutput(ns("dea_stat_para"))
                  ),
                  br(),
                  div(
                    style = "",
                    tableOutput(ns("DEA_genes"))
                  ),
                  br(),
                  downloadButton(ns("DEA_genes_export"), "Download full table", class = "btn btn-primary btn-sm")
                )
              ),
              
              tabPanel(
                "Feature Plots of Top Differentially Expressed Genes",
                
                box(
                  width = NULL,
                  height = NULL,
                  solidHeader = TRUE,
                  status = "primary",
                  uiOutput(ns("group_choice")),
                  br(),
                  div(
                    style = "padding-left: 6em; font-size: 12pt; color: black;",
                    textOutput(ns("fea_plot_stat"))
                  ),
                  br(),
                  plotlyOutput(ns("TopDEA"), height = 600),
                  br(),
                  radioButtons(
                    ns("fmt_fea_plot"),
                    "Download format",
                    choices = c("PDF" = "pdf", "SVG" = "svg"),
                    inline = TRUE
                  ),
                  downloadButton(ns("fea_plot_export"), "Download feature plot", class = "btn btn-primary btn-sm")
                )
              ),
              
              tabPanel(
                "GO Biological Process",
                box(
                  width = NULL,
                  height = NULL,
                  solidHeader = TRUE,
                  status = "primary",
                  plotlyOutput(ns("GO_BP"), height = 600),
                  br(),
                  radioButtons(
                    ns("fmt_GO_BP"),
                    "Download format",
                    choices = c("PDF" = "pdf", "SVG" = "svg"),
                    inline = TRUE
                  ),
                  downloadButton(ns("fea_GO_BP_export"), "Download GO BP plot", class = "btn btn-primary btn-sm"),
                  downloadButton(ns("download_GO_BP"), "Download GO BP table", class = "btn btn-primary btn-sm")
                )
              ),
              
              tabPanel(
                "GO Molecular Function",
                box(
                  width = NULL,
                  height = NULL,
                  solidHeader = TRUE,
                  status = "primary",
                  plotlyOutput(ns("GO_MF"), height = 600),
                  br(),
                  radioButtons(
                    ns("fmt_GO_MF"),
                    "Download format",
                    choices = c("PDF" = "pdf", "SVG" = "svg"),
                    inline = TRUE
                  ),
                  downloadButton(ns("fea_GO_MF_export"), "Download GO MF plot", class = "btn btn-primary btn-sm"),
                  downloadButton(ns("download_GO_MF"), "Download GO MF table", class = "btn btn-primary btn-sm")
                )
              ),
              
              tabPanel(
                "GO Cellular Component",
                box(
                  width = NULL,
                  height = NULL,
                  solidHeader = TRUE,
                  status = "primary",
                  plotlyOutput(ns("GO_CC"), height = 600),
                  br(),
                  radioButtons(
                    ns("fmt_GO_CC"),
                    "Download format",
                    choices = c("PDF" = "pdf", "SVG" = "svg"),
                    inline = TRUE
                  ),
                  downloadButton(ns("fea_GO_CC_export"), "Download GO CC plot", class = "btn btn-primary btn-sm"),
                  downloadButton(ns("download_GO_CC"), "Download GO CC table", class = "btn btn-primary btn-sm")
                )
              ),
              
              tabPanel(
                "KEGG",
                box(
                  width = NULL,
                  height = NULL,
                  solidHeader = TRUE,
                  status = "primary",
                  plotlyOutput(ns("KEGG"), height = 600),
                  br(),
                  radioButtons(
                    ns("fmt_KEGG"),
                    "Download format",
                    choices = c("PDF" = "pdf", "SVG" = "svg"),
                    inline = TRUE
                  ),
                  downloadButton(ns("fea_KEGG_export"), "Download KEGG plot", class = "btn btn-primary btn-sm"),
                  downloadButton(ns("download_KEGG"), "Download KEGG table", class = "btn btn-primary btn-sm")
                )
              )
            )
          )
        )
      )
    )
    
    
  )
}


#------------------------------------------
#    DE Expression server
#-----------------------------------------

mod_sc_de_server <- function(id, rval) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    params <- reactiveValues(
      df_dea = NULL,
      catch_param_vis = NULL,
      catch_param_dea = NULL,
      check = TRUE,
      comment_out = NULL,
      OrgDb = "org.Hs.eg.db",
      organism = "hsa",
      ego_BP = NULL,
      ego_MF = NULL,
      ego_CC = NULL,
      KEGG = NULL
    )
    
    observeEvent(input[["HELP_DEG_INPUT"]], {
      showModal(modalDialog(
        title = HTML(help_infos[help_infos$key == "HELP_DEG_INPUT", "title"]),
        HTML(help_infos[help_infos$key == "HELP_DEG_INPUT", "value"]),
        easyClose = TRUE,
        footer = modalButton("OK")
      ))
    })
    
    observeEvent(input[["HELP_DEG_PARAMS"]], {
      showModal(modalDialog(
        title = HTML(help_infos[help_infos$key == "HELP_DEG_PARAMS", "title"]),
        HTML(help_infos[help_infos$key == "HELP_DEG_PARAMS", "value"]),
        easyClose = TRUE,
        footer = modalButton("OK")
      ))
    })
    
    observeEvent(input[["HELP_DEG_VIS"]], {
      showModal(modalDialog(
        title = HTML(help_infos[help_infos$key == "HELP_DEG_VIS", "title"]),
        HTML(help_infos[help_infos$key == "HELP_DEG_VIS", "value"]),
        easyClose = TRUE,
        footer = modalButton("OK")
      ))
    })
    
    observeEvent(input[["HELP_DEG_TABLE"]], {
      showModal(modalDialog(
        title = HTML(help_infos[help_infos$key == "HELP_DEG_TABLE", "title"]),
        HTML(help_infos[help_infos$key == "HELP_DEG_TABLE", "value"]),
        easyClose = TRUE,
        footer = modalButton("OK")
      ))
    })
    
    numericCol <- reactive({
      if (is.null(rval$seurat)) return(NULL)
      colnames(dplyr::select_if(rval$seurat@meta.data, is.numeric))
    })
    
    factorCol <- reactive({
      if (is.null(rval$seurat)) return(NULL)
      colnames(dplyr::select_if(rval$seurat@meta.data, is.factor))
    })
    
    getGenes <- reactive({
      if (is.null(rval$seurat)) return(NULL)
      rownames(rval$seurat)
    })
    
    output$select_meta_fac1 <- renderUI({
      if (is.null(rval$seurat)) return(NULL)
      selectizeInput(
        ns("selectMetafac1"),
        "Choose Annotation for Group 1",
        factorCol(),
        selected = NULL,
        multiple = TRUE,
        options = list(maxItems = 1, placeholder = "select 1 metadata max")
      )
    })
    
    output$range_meta_fac1 <- renderUI({
      req(input$selectMetafac1)
      selectizeInput(
        ns("selectLevelsFac1"),
        "Keep Labels for Group 1",
        choices = levels(droplevels(rval$seurat@meta.data[, input$selectMetafac1[1]])),
        selected = levels(droplevels(rval$seurat@meta.data[, input$selectMetafac1[1]])),
        multiple = TRUE
      )
    })
    
    output$select_meta_fac2 <- renderUI({
      if (is.null(rval$seurat) ||
          is.null(input$selectMetafac1) ||
          is.null(input$selectLevelsFac1)) {
        return(NULL)
      }
      selectizeInput(
        ns("selectMetafac2"),
        "Choose Annotation for Group 2",
        factorCol(),
        selected = NULL,
        multiple = TRUE,
        options = list(maxItems = 1, placeholder = "select 1 metadata max")
      )
    })
    
    output$range_meta_fac2 <- renderUI({
      req(input$selectMetafac2)
      if (length(input$selectMetafac2) != 1) return(NULL)
      if (is.null(input$selectLevelsFac1)) return(NULL)
      selectizeInput(
        ns("selectLevelsFac2"),
        "Keep Labels for Group 2",
        choices = levels(droplevels(rval$seurat@meta.data[, input$selectMetafac2[1]])),
        selected = levels(droplevels(rval$seurat@meta.data[, input$selectMetafac2[1]])),
        multiple = TRUE
      )
    })
    
    output$select_reduc <- renderUI({
      if (is.null(input$selectMetafac1) || is.null(input$selectMetafac2)) return(NULL)
      reduc <- Seurat::Reductions(rval$seurat)
      if (length(reduc) == 0) return(NULL)
      selectInput(
        ns("sel_rd_vis"),
        label = "Choose Embedding for Display",
        choices = reduc,
        selected = tail(reduc, 1)
      )
    })
    
    output$dea_vis_btn <- renderUI({
      if (is.null(input$selectMetafac1) && is.null(input$selectMetafac2)) return(NULL)
      actionButton(ns("dea_vis_btn2"), label = "Next")
    })
    
    output$select_method_Markers <- renderUI({
      if (is.null(rval$seurat)) return(NULL)
      selectInput(
        ns("method_Markers"),
        label = "Choose Method to Identify Differentially Expressed Genes",
        choices = list(
          "wilcox" = "wilcox",
          "bimod" = "bimod",
          "t" = "t",
          "negbinom" = "negbinom",
          "poisson" = "poisson",
          "LR" = "LR"
        ),
        selected = "wilcox"
      )
    })
    
    output$min_pct <- renderUI({
      numericInput(
        ns("min"),
        "Set the Minimal Fraction of Cells Expressing Genes",
        value = 0.25,
        min = 0,
        max = 1,
        step = 0.1
      )
    })
    
    output$nb_logfc <- renderUI({
      numericInput(
        ns("logfc"),
        "Set the Minimal logFC of Gene Expression Difference",
        value = 0.25,
        min = 0,
        step = 0.1
      )
    })
    
    output$checkbox_only_pos <- renderUI({
      checkboxInput(
        ns("only_pos"),
        "Select only positive markers",
        value = TRUE
      )
    })
    
    output$dea_run_btn <- renderUI({
      actionButton(ns("dea_run_btn2"), label = "Run Analysis")
    })
    
    dea_vis_plot_gg <- eventReactive(input$dea_vis_btn2, {
      req(input$sel_rd_vis, input$selectMetafac1, input$selectMetafac2)
      
      reduc <- input$sel_rd_vis
      
      fac_sel1 <- input$selectMetafac1
      fac_sel2 <- input$selectMetafac2
      
      para_sel1 <- input$selectLevelsFac1
      para_sel2 <- input$selectLevelsFac2
      
      group1 <- which(rval$seurat@meta.data[, fac_sel1] %in% para_sel1)
      group2 <- which(rval$seurat@meta.data[, fac_sel2] %in% para_sel2)
      common_cells <- intersect(group1, group2)
      
      rval$seurat@meta.data$selected_groups <- "Others"
      rval$seurat@meta.data$selected_groups[group1] <- "Group1"
      rval$seurat@meta.data$selected_groups[group2] <- "Group2"
      
      if (length(common_cells) > 0) {
        rval$seurat@meta.data$selected_groups[common_cells] <- "Common"
        
        shinyalert(
          "Warning",
          paste(
            "Groups with common cells have been selected. Common cells (N) =",
            length(common_cells)
          ),
          type = "warning"
        )
        
        groups <- c("Group1", "Group2")
        sel_groups <- unique(rval$seurat@meta.data[, "selected_groups"])
        
        if (length(intersect(groups, sel_groups)) < 2) {
          shinyalert(
            "Warning",
            "One of the groups is a subset of the other group. Please choose other groups",
            type = "warning"
          )
        }
      }
      
      para_st_1 <- paste(para_sel1, collapse = "-")
      para_st_2 <- paste(para_sel2, collapse = "-")
      
      params$catch_param_vis <- paste(
        "Group 1 (", fac_sel1, " [", para_st_1, "]) and Group 2 (",
        fac_sel2, " [", para_st_2, "])",
        sep = ""
      )
      
      if (length(common_cells) > 0) {
        output$state_vis_out <- renderText({
          paste(
            "Common cells (will be excluded) (N):", length(common_cells),
            "; Group 1 (", paste(para_sel1, collapse = ","), "):", length(group1),
            "; Group 2 (", paste(para_sel2, collapse = ","), "):", length(group2)
          )
        })
      } else {
        output$state_vis_out <- renderText({
          paste(
            "Group 1 (", paste(para_sel1, collapse = ","), "):", length(group1),
            "; Group 2 (", paste(para_sel2, collapse = ","), "):", length(group2)
          )
        })
      }
      
      params$check <- FALSE
      
      Seurat::DimPlot(
        rval$seurat,
        reduction = reduc,
        dims = c(1, 2),
        group.by = "selected_groups",
        pt.size = 0.5,
        label = FALSE
      )
    })
    
    output$dea_vis_plot <- renderPlotly({
      req(dea_vis_plot_gg())
      plotly::ggplotly(dea_vis_plot_gg())
    })
    
    output$down_dea_vis <- downloadHandler(
      filename = function() {
        paste0("DEG_defined_groups.", input$fmt_dea_vis)
      },
      content = function(file) {
        req(dea_vis_plot_gg())
        if (input$fmt_dea_vis == "pdf") {
          grDevices::pdf(file, width = 7, height = 7)
          print(dea_vis_plot_gg())
          grDevices::dev.off()
        } else if (input$fmt_dea_vis == "svg") {
          grDevices::svg(file, width = 7, height = 7)
          print(dea_vis_plot_gg())
          grDevices::dev.off()
        }
      }
    )
    
    observeEvent(input$dea_run_btn2, {
      updateTabsetPanel(
        session = session,
        inputId = "panel_deg",
        selected = "Table of Differentially Expressed Genes"
      )
    })
    
    DEA_table <- eventReactive(input$dea_run_btn2, {
      params$df_dea <- NULL
      
      if (!is.null(rval$seurat@meta.data[["species"]]) &&
          rval$seurat@meta.data[["species"]][1] == "mouse") {
        params$OrgDb <- "org.Mm.eg.db"
        params$organism <- "mmu"
      }
      
      fac_sel1 <- input$selectMetafac1
      fac_sel2 <- input$selectMetafac2
      
      para_sel1 <- input$selectLevelsFac1
      para_sel2 <- input$selectLevelsFac2
      
      Seurat::Idents(rval$seurat) <- rval$seurat@meta.data$selected_groups
      
      if ("Group1" %in% levels(Seurat::Idents(rval$seurat)) &&
          "Group2" %in% levels(Seurat::Idents(rval$seurat)) &&
          table(Seurat::Idents(rval$seurat))["Group1"] > 3 &&
          table(Seurat::Idents(rval$seurat))["Group2"] > 3) {
        
        show_modal_spinner(text = "Markers finding. Please wait", spin = "circle")
        
        temp_df_dea <- Seurat::FindMarkers(
          object = rval$seurat,
          ident.1 = "Group1",
          ident.2 = "Group2",
          group.by = "selected_groups",
          test.use = input$method_Markers,
          only.pos = input$only_pos,
          min.pct = input$min,
          logfc.threshold = input$logfc
        )
        
        remove_modal_spinner()
      } else {
        showModal(modalDialog(
          "Re run the differential expression analysis for newly defined groups",
          size = "s",
          easyClose = TRUE,
          footer = modalButton("OK")
        ))
        return(NULL)
      }
      
      if (nrow(temp_df_dea) <= 0) {
        showModal(modalDialog(
          "No features pass logfc threshold",
          size = "s",
          easyClose = TRUE,
          footer = modalButton("OK")
        ))
        return(NULL)
      }
      
      params$catch_param_dea <- paste(
        "Chosen Statistical Method:", input$method_Markers,
        "; Only positive genes:", input$only_pos,
        "; Minimal Fraction of Cells Expressing Genes:", input$min,
        "; Minimal logFC:", input$logfc
      )
      
      temp_df_dea$diff.pct <- temp_df_dea$pct.1 - temp_df_dea$pct.2
      temp_df_dea <- temp_df_dea[, c("p_val", "avg_log2FC", "pct.1", "pct.2", "diff.pct", "p_val_adj")]
      temp_df_dea$p_val <- as.numeric(scales::scientific(temp_df_dea$p_val, digits = 3))
      temp_df_dea$p_val_adj <- as.numeric(scales::scientific(temp_df_dea$p_val_adj, digits = 3))
      temp_df_dea[, c("avg_log2FC", "pct.1", "pct.2", "diff.pct")] <-
        round(temp_df_dea[, c("avg_log2FC", "pct.1", "pct.2", "diff.pct")], 3)
      
      params$check <- TRUE
      params$df_dea <- temp_df_dea
      
      temp_df_dea
    })
    
    output$DEA_genes <- function() {
      if (is.null(DEA_table())) return(NULL)
      if (nrow(DEA_table()) > 200) {
        knitr::kable(
          DEA_table()[1:200, ],
          "html",
          caption = "Only top 200 genes are shown. Download to see all genes"
        ) |>
          kableExtra::kable_styling(bootstrap_options = c("striped", "hover")) |>
          kableExtra::scroll_box(width = "100%", height = "450px")
      } else {
        knitr::kable(DEA_table(), "html") |>
          kableExtra::kable_styling(bootstrap_options = c("striped", "hover")) |>
          kableExtra::scroll_box(width = "100%", height = "450px")
      }
    }
    
    output$DEA_genes_export <- downloadHandler(
      filename = function() {
        paste0("DEA_genes_", Sys.Date(), ".csv")
      },
      content = function(file) {
        req(DEA_table())
        write.csv(DEA_table(), file, row.names = TRUE)
      }
    )
    
    observeEvent(input$dea_run_btn2, {
      output$dea_table_stat <- renderUI({
        if (params$check) {
          params$comment_out <- paste(
            "Differentially expressed genes for",
            params$catch_param_vis,
            ";",
            params$catch_param_dea
          )
          str1 <- paste("Differentially expressed genes for", params$catch_param_vis)
          str2 <- params$catch_param_dea
          HTML(paste(str1, str2, sep = "<br/><br/>"))
        } else {
          "Re run the differential expression analysis for newly defined groups"
        }
      })
    })
    
    DEA_featurePlot_gg <- eventReactive(params$df_dea, {
      req(params$df_dea, input$sel_rd_vis)
      df <- params$df_dea
      df <- df[order(df$p_val_adj, abs(df$avg_log2FC), decreasing = c(FALSE, TRUE)), ]
      genes <- rownames(df)
      dims <- c(1, 2)
      
      if (length(genes) > 9) {
        Seurat::FeaturePlot(
          rval$seurat,
          reduction = input$sel_rd_vis,
          dims = dims,
          features = genes[1:9],
          ncol = 3,
          pt.size = 0.5
        )
      } else if (length(genes) == 0) {
        NULL
      } else {
        Seurat::FeaturePlot(
          rval$seurat,
          reduction = input$sel_rd_vis,
          dims = dims,
          features = genes,
          ncol = 3,
          pt.size = 0.5
        )
      }
    })
    
    output$TopDEA <- renderPlotly({
      req(DEA_featurePlot_gg())
      plotly::ggplotly(DEA_featurePlot_gg())
    })
    
    output$fea_plot_export <- downloadHandler(
      filename = function() {
        paste0("top_differentially_expressed_genes.", input$fmt_fea_plot)
      },
      content = function(file) {
        req(DEA_featurePlot_gg())
        if (input$fmt_fea_plot == "pdf") {
          grDevices::pdf(file, width = 7, height = 7)
          print(DEA_featurePlot_gg())
          grDevices::dev.off()
        } else if (input$fmt_fea_plot == "svg") {
          grDevices::svg(file, width = 7, height = 7)
          print(DEA_featurePlot_gg())
          grDevices::dev.off()
        }
      }
    )
    
    observeEvent(input$dea_run_btn2, {
      output$fea_plot_stat <- renderText({
        if (params$check) {
          paste(
            "Differentially expressed genes for",
            params$catch_param_vis,
            ";",
            params$catch_param_dea
          )
        } else {
          "Re run the differential expression analysis for newly defined groups"
        }
      })
    })
    
    GO_BP_Plot_gg <- eventReactive(params$df_dea, {
      req(params$df_dea)
      df <- params$df_dea
      df <- df[order(df$p_val_adj, abs(df$avg_log2FC), decreasing = c(FALSE, TRUE)), ]
      genes <- rownames(df)
      
      eg <- clusterProfiler::bitr(
        genes,
        fromType = "SYMBOL",
        toType = "ENTREZID",
        OrgDb = params$OrgDb,
        drop = FALSE
      )
      eg <- eg[!duplicated(eg["SYMBOL"]), ]
      eg <- stats::na.omit(eg)
      
      if (nrow(eg) == 0) {
        params$ego_BP <- NULL
        showModal(modalDialog(
          "No significant result",
          size = "s",
          easyClose = TRUE,
          footer = modalButton("OK")
        ))
        return(NULL)
      }
      
      show_modal_spinner(text = "Enrichment processing. Please wait", spin = "circle")
      
      ego_BP <- clusterProfiler::enrichGO(
        gene = eg[, 2],
        OrgDb = params$OrgDb,
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.05,
        readable = TRUE
      )
      
      remove_modal_spinner()
      
      if (!is.null(ego_BP) && nrow(ego_BP@result) > 0) {
        params$ego_BP <- ego_BP@result
        enrichplot::dotplot(ego_BP, showCategory = 30)
      } else {
        params$ego_BP <- NULL
        showModal(modalDialog(
          "No significant result",
          size = "s",
          easyClose = TRUE,
          footer = modalButton("OK")
        ))
        NULL
      }
    })
    
    output$GO_BP <- renderPlotly({
      req(GO_BP_Plot_gg())
      plotly::ggplotly(GO_BP_Plot_gg())
    })
    
    output$fea_GO_BP_export <- downloadHandler(
      filename = function() {
        paste0("GO_Biological_Process_plot.", input$fmt_GO_BP)
      },
      content = function(file) {
        req(GO_BP_Plot_gg())
        if (input$fmt_GO_BP == "pdf") {
          grDevices::pdf(file, width = 7, height = 7)
          print(GO_BP_Plot_gg())
          grDevices::dev.off()
        } else if (input$fmt_GO_BP == "svg") {
          grDevices::svg(file, width = 7, height = 7)
          print(GO_BP_Plot_gg())
          grDevices::dev.off()
        }
      }
    )
    
    output$download_GO_BP <- downloadHandler(
      filename = function() {
        paste("GO_Biological_Process_Table", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        con <- file(file, open = "wt")
        writeLines("# Results were obtained from CoTRA", con)
        if (!is.null(params$ego_BP)) {
          write.csv(params$ego_BP, con, row.names = FALSE)
        }
        close(con)
      }
    )
    
    GO_MF_Plot_gg <- eventReactive(params$df_dea, {
      req(params$df_dea)
      df <- params$df_dea
      df <- df[order(df$p_val_adj, abs(df$avg_log2FC), decreasing = c(FALSE, TRUE)), ]
      genes <- rownames(df)
      
      eg <- clusterProfiler::bitr(
        genes,
        fromType = "SYMBOL",
        toType = "ENTREZID",
        OrgDb = params$OrgDb,
        drop = FALSE
      )
      eg <- eg[!duplicated(eg["SYMBOL"]), ]
      eg <- stats::na.omit(eg)
      
      if (nrow(eg) == 0) {
        params$ego_MF <- NULL
        showModal(modalDialog(
          "No significant result",
          size = "s",
          easyClose = TRUE,
          footer = modalButton("OK")
        ))
        return(NULL)
      }
      
      show_modal_spinner(text = "Enrichment processing. Please wait", spin = "circle")
      
      ego_MF <- clusterProfiler::enrichGO(
        gene = eg[, 2],
        OrgDb = params$OrgDb,
        ont = "MF",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.05,
        readable = TRUE
      )
      
      remove_modal_spinner()
      
      if (!is.null(ego_MF) && nrow(ego_MF@result) > 0) {
        params$ego_MF <- ego_MF@result
        enrichplot::dotplot(ego_MF, showCategory = 30)
      } else {
        params$ego_MF <- NULL
        showModal(modalDialog(
          "No significant result",
          size = "s",
          easyClose = TRUE,
          footer = modalButton("OK")
        ))
        NULL
      }
    })
    
    output$GO_MF <- renderPlotly({
      req(GO_MF_Plot_gg())
      plotly::ggplotly(GO_MF_Plot_gg())
    })
    
    output$fea_GO_MF_export <- downloadHandler(
      filename = function() {
        paste0("GO_Molecular_Function_plot.", input$fmt_GO_MF)
      },
      content = function(file) {
        req(GO_MF_Plot_gg())
        if (input$fmt_GO_MF == "pdf") {
          grDevices::pdf(file, width = 7, height = 7)
          print(GO_MF_Plot_gg())
          grDevices::dev.off()
        } else if (input$fmt_GO_MF == "svg") {
          grDevices::svg(file, width = 7, height = 7)
          print(GO_MF_Plot_gg())
          grDevices::dev.off()
        }
      }
    )
    
    output$download_GO_MF <- downloadHandler(
      filename = function() {
        paste("GO_Molecular_Function_Table", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        con <- file(file, open = "wt")
        writeLines("# Results were obtained from CoTRA", con)
        if (!is.null(params$ego_MF)) {
          write.csv(params$ego_MF, con, row.names = FALSE)
        }
        close(con)
      }
    )
    
    GO_CC_Plot_gg <- eventReactive(params$df_dea, {
      req(params$df_dea)
      df <- params$df_dea
      df <- df[order(df$p_val_adj, abs(df$avg_log2FC), decreasing = c(FALSE, TRUE)), ]
      genes <- rownames(df)
      
      eg <- clusterProfiler::bitr(
        genes,
        fromType = "SYMBOL",
        toType = "ENTREZID",
        OrgDb = params$OrgDb,
        drop = FALSE
      )
      eg <- eg[!duplicated(eg["SYMBOL"]), ]
      eg <- stats::na.omit(eg)
      
      if (nrow(eg) == 0) {
        params$ego_CC <- NULL
        showModal(modalDialog(
          "No significant result",
          size = "s",
          easyClose = TRUE,
          footer = modalButton("OK")
        ))
        return(NULL)
      }
      
      show_modal_spinner(text = "Enrichment processing. Please wait", spin = "circle")
      
      ego_CC <- clusterProfiler::enrichGO(
        gene = eg[, 2],
        OrgDb = params$OrgDb,
        ont = "CC",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.05,
        readable = TRUE
      )
      
      remove_modal_spinner()
      
      if (!is.null(ego_CC) && nrow(ego_CC@result) > 0) {
        params$ego_CC <- ego_CC@result
        enrichplot::dotplot(ego_CC, showCategory = 30)
      } else {
        params$ego_CC <- NULL
        showModal(modalDialog(
          "No significant result",
          size = "s",
          easyClose = TRUE,
          footer = modalButton("OK")
        ))
        NULL
      }
    })
    
    output$GO_CC <- renderPlotly({
      req(GO_CC_Plot_gg())
      plotly::ggplotly(GO_CC_Plot_gg())
    })
    
    output$fea_GO_CC_export <- downloadHandler(
      filename = function() {
        paste0("GO_Cellular_Component_plot.", input$fmt_GO_CC)
      },
      content = function(file) {
        req(GO_CC_Plot_gg())
        if (input$fmt_GO_CC == "pdf") {
          grDevices::pdf(file, width = 7, height = 7)
          print(GO_CC_Plot_gg())
          grDevices::dev.off()
        } else if (input$fmt_GO_CC == "svg") {
          grDevices::svg(file, width = 7, height = 7)
          print(GO_CC_Plot_gg())
          grDevices::dev.off()
        }
      }
    )
    
    output$download_GO_CC <- downloadHandler(
      filename = function() {
        paste("GO_Cellular_Component_Table", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        con <- file(file, open = "wt")
        writeLines("# Results were obtained from CoTRA", con)
        if (!is.null(params$ego_CC)) {
          write.csv(params$ego_CC, con, row.names = FALSE)
        }
        close(con)
      }
    )
    
    KEGG_Plot_gg <- eventReactive(params$df_dea, {
      req(params$df_dea)
      df <- params$df_dea
      df <- df[order(df$p_val_adj, abs(df$avg_log2FC), decreasing = c(FALSE, TRUE)), ]
      genes <- rownames(df)
      
      eg <- clusterProfiler::bitr(
        genes,
        fromType = "SYMBOL",
        toType = "ENTREZID",
        OrgDb = params$OrgDb,
        drop = FALSE
      )
      eg <- eg[!duplicated(eg["SYMBOL"]), ]
      eg <- stats::na.omit(eg)
      
      if (nrow(eg) == 0) {
        params$KEGG <- NULL
        showModal(modalDialog(
          "No significant result",
          size = "s",
          easyClose = TRUE,
          footer = modalButton("OK")
        ))
        return(NULL)
      }
      
      show_modal_spinner(text = "Enrichment processing. Please wait", spin = "circle")
      
      KEGG <- clusterProfiler::enrichKEGG(
        gene = eg[, 2],
        organism = params$organism,
        pvalueCutoff = 0.05
      )
      
      remove_modal_spinner()
      
      if (!is.null(KEGG) && nrow(KEGG@result) > 0) {
        params$KEGG <- KEGG@result
        enrichplot::dotplot(KEGG, showCategory = 30)
      } else {
        params$KEGG <- NULL
        showModal(modalDialog(
          "No significant result",
          size = "s",
          easyClose = TRUE,
          footer = modalButton("OK")
        ))
        NULL
      }
    })
    
    output$KEGG <- renderPlotly({
      req(KEGG_Plot_gg())
      plotly::ggplotly(KEGG_Plot_gg())
    })
    
    output$fea_KEGG_export <- downloadHandler(
      filename = function() {
        paste0("KEGG_plot.", input$fmt_KEGG)
      },
      content = function(file) {
        req(KEGG_Plot_gg())
        if (input$fmt_KEGG == "pdf") {
          grDevices::pdf(file, width = 7, height = 7)
          print(KEGG_Plot_gg())
          grDevices::dev.off()
        } else if (input$fmt_KEGG == "svg") {
          grDevices::svg(file, width = 7, height = 7)
          print(KEGG_Plot_gg())
          grDevices::dev.off()
        }
      }
    )
    
    output$download_KEGG <- downloadHandler(
      filename = function() {
        paste("KEGG_Table", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        con <- file(file, open = "wt")
        writeLines("# Results were obtained from CoTRA", con)
        if (!is.null(params$KEGG)) {
          write.csv(params$KEGG, con, row.names = FALSE)
        }
        close(con)
      }
    )
    
    rval
    
    
  })
}
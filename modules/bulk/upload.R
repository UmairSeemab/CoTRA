mod_bulk_upload_ui <- function(id){
  ns <- NS(id)
  
  tagList(
    shinyjs::useShinyjs(),
    
    tags$script('$(function () {$("[data-toggle=\'tooltip\']").tooltip();});'),
    
    actionButton(ns("help_toggle"), icon("question-circle"), class = "btn btn-info btn-sm"),
    
    hidden(
      div(
        id = ns("help_box"),
        style = "border:1px solid #d9d9d9; padding:12px; margin-top:10px; border-radius:6px; background:#f7f7f7;",
        
        h4("Step-by-Step Guide"),
        h5("1. Upload counts"),
        p("Upload a raw counts file (by clicking browse button and selecting file from local directory) separated by comma, semicolon, or tab."),
        h5("2. Confirm detected gene ID column"),
        p("CoTRA identifies the most likely ID column. Check it."),
        h5("3. Choose normalization method"),
        p("Choose FPKM, RPKM, TPM or TMM depending on your analysis plan."),
        h5("4. Load data"),
        p("Press load data button to load the gene count file."),
        h5("5. Confirm gene id column"),
        p("After confirming, you see sample previews, gene previews, and group suggestions.")
      )
    ),
    
    br(),
    
    radioGroupButtons(
      ns("data_type"),
      "Input data type",
      choices = c("Counts matrix" = "counts"),
      justified = TRUE
    ),
    
    fileInput(
      ns("matrix_file"),
      "Upload matrix file",
      accept = c(".txt", ".tsv", ".csv", ".tab")
    ),
    
    selectInput(
      ns("delim_choice"),
      "Separator",
      choices = c("Auto detect" = "auto",
                  "Comma" = "comma",
                  "Semicolon" = "semicolon",
                  "Tab" = "tab"),
      selected = "auto"
    ),
    
    textInput(ns("id_col"), "Gene ID column (optional)", value = ""),
    
    radioGroupButtons(
      ns("norm_method"),
      "Normalization",
      choices = c("None" = "none",
                  "FPKM" = "fpkm",
                  "RPKM" = "rpkm",
                  "TPM" = "tpm",
                  "TMM" = "tmm"),
      selected = "none",
      justified = TRUE,
      status = "primary",
      checkIcon = list(yes = icon("check"))
    ),
    
    div(
      actionButton(ns("load_btn"), "Load data", class = "btn btn-primary"),
      HTML(paste0(
        '<i id="', ns("tt_load"),
        '" class="fa fa-question-circle"
        style="margin-left:6px; cursor:pointer;"
        data-toggle="tooltip" data-placement="right"
        title="Loads your count matrix and detects gene ID column"></i>'
      ))
    ),
    
    br(), br(),
    
    uiOutput(ns("id_detect_ui")),
    
    verbatimTextOutput(ns("summary")),
    strong("Sample names"),
    tableOutput(ns("sample_preview")),
    
    br(),
    strong("Gene IDs preview"),
    tableOutput(ns("gene_preview")),
    
    br(),
    strong("Suggested groups"),
    verbatimTextOutput(ns("group_suggest")),
    
    uiOutput(ns("download_ui"))
  )
}


mod_bulk_upload_server <- function(id){
  moduleServer(id, function(input, output, session){
    
    ns <- session$ns
    library(shinyjs)
    library(readr)
    library(edgeR)
    library(DESeq2)
    
    shinyjs::useShinyjs()
    
    observeEvent(input$help_toggle, {
      shinyjs::toggle(id = "help_box", anim = TRUE)
    })
    
    
    detect_gene_column <- function(df){
      cn <- colnames(df)
      if (length(cn) == 0) return(1)
      if (cn[1] == "" || grepl("^\\.\\.\\.", cn[1])) return(1)
      
      score <- sapply(df, function(col){
        vals <- as.character(col[1:min(20, length(col))])
        ens <- mean(grepl("^ENS", vals))
        symbol <- mean(grepl("^[A-Za-z0-9\\-]+$", vals))
        ens * 2 + symbol
      })
      
      which.max(score)
    }
    
    infer_delim <- function(file){
      first_line <- readLines(file, n = 1, warn = FALSE)
      if (!length(first_line)) return("tab")
      
      counts <- c(
        comma = lengths(regmatches(first_line, gregexpr(",", first_line))),
        semicolon = lengths(regmatches(first_line, gregexpr(";", first_line))),
        tab = lengths(regmatches(first_line, gregexpr("\t", first_line)))
      )
      
      best <- names(which.max(counts))
      if (is.null(best) || is.na(best) || counts[[best]] == 0) return("tab")
      best
    }
    
    delim_to_char <- function(key){
      if (key == "comma") return(",")
      if (key == "semicolon") return(";")
      "\t"
    }
    
    safe_numeric <- function(mat){
      suppressWarnings(apply(mat, 2, function(x) as.numeric(x)))
    }
    
    guess_groups <- function(samples){
      s <- tolower(samples)
      case <- grepl("tumor|treated|case|ko|exp|male", s)
      ctrl <- grepl("control|wt|normal|untreated|vehicle|female", s)
      if (any(case) | any(ctrl)) {
        return(paste0(
          "Case: ", paste(samples[case], collapse = ", "),
          "\nControl: ", paste(samples[ctrl], collapse = ", ")
        ))
      }
      "No pattern detected"
    }
    
    normalize_counts <- function(counts, method){
      if (method == "none") return(counts)
      
      gene_length <- rep(1000, nrow(counts))
      
      if (method %in% c("rpkm", "fpkm")) {
        return(edgeR::rpkm(counts, gene.length = gene_length))
      }
      
      if (method == "tpm") {
        rpk <- counts / (gene_length / 1000)
        return(t(t(rpk) / colSums(rpk) * 1e6))
      }
      
      if (method == "tmm") {
        y <- DGEList(counts = counts)
        y <- calcNormFactors(y)
        return(cpm(y, log = FALSE))
      }
      
      counts
    }
    
    
    rv <- reactiveValues(
      raw = NULL,
      norm = NULL,
      final = NULL,
      df = NULL,
      detected_col = NULL
    )
    
    
    observeEvent(input$load_btn, {
      req(input$matrix_file)
      
      detected <- infer_delim(input$matrix_file$datapath)
      chosen <- input$delim_choice
      if (is.null(chosen) || chosen == "auto") {
        delim_key <- detected
      } else {
        delim_key <- chosen
      }
      
      df <- read_delim(
        input$matrix_file$datapath,
        delim = delim_to_char(delim_key),
        show_col_types = FALSE,
        trim_ws = TRUE,
        progress = FALSE
      )
      
      df <- df[, colSums(!is.na(df)) > 0, drop = FALSE]
      
      detected_col_index <- detect_gene_column(df)
      detected_col_name  <- colnames(df)[detected_col_index]
      
      rv$detected_col <- detected_col_index
      rv$df <- df
      
      output$id_detect_ui <- renderUI({
        div(
          h4("Detected Gene ID column"),
          p("Detected:", strong(detected_col_name)),
          p("Separator used:", strong(delim_key)),
          
          tableOutput(ns("col_preview")),
          
          div(
            actionButton(ns("confirm_id"), "Confirm Gene ID column", class = "btn btn-success"),
            HTML(paste0(
              '<i id="', ns("tt_confirm"),
              '" class="fa fa-question-circle"
              style="margin-left:6px; cursor:pointer;"
              data-toggle="tooltip" data-placement="right"
              title="Locks in the detected ID column and converts matrix to numeric values"></i>'
            ))
          )
        )
      })
      
      output$col_preview <- renderTable({
        head(df[[detected_col_name]], 10)
      })
    })
    
    
    observeEvent(input$confirm_id, {
      req(rv$df)
      req(rv$detected_col)
      
      df <- rv$df
      
      manual_id <- trimws(input$id_col)
      if (nzchar(manual_id) && manual_id %in% colnames(df)) {
        id_col <- manual_id
      } else {
        id_col <- colnames(df)[rv$detected_col]
      }
      
      gene_ids <- df[[id_col]]
      df[[id_col]] <- NULL
      
      m <- safe_numeric(as.matrix(df))
      rownames(m) <- gene_ids
      
      rv$raw <- m
      
      norm_m <- normalize_counts(m, input$norm_method)
      rv$norm <- norm_m
      rv$final <- norm_m
      
      output$id_detect_ui <- renderUI(NULL)
      
      output$summary <- renderText({
        paste0(
          "Genes: ", nrow(norm_m),
          "\nSamples: ", ncol(norm_m),
          "\nNormalization: ", input$norm_method
        )
      })
      
      output$sample_preview <- renderTable({
        head(data.frame(samples = colnames(norm_m)), 10)
      })
      
      output$gene_preview <- renderTable({
        head(data.frame(GeneID = rownames(norm_m)), 10)
      })
      
      output$group_suggest <- renderText({
        guess_groups(colnames(norm_m))
      })
      
      output$download_ui <- renderUI({
        if (input$norm_method != "none") {
          div(
            downloadButton(ns("download_norm"), "Download normalized CSV"),
            HTML(paste0(
              '<i id="', ns("tt_download"),
              '" class="fa fa-question-circle"
              style="margin-left:6px; cursor:pointer;"
              data-toggle="tooltip" data-placement="right"
              title="Saves the normalized count matrix as a CSV file"></i>'
            ))
          )
        } else {
          NULL
        }
      })
      
      output$download_norm <- downloadHandler(
        filename = function(){ "normalized_matrix.csv" },
        content = function(file){
          write.csv(rv$norm, file, row.names = TRUE)
        }
      )
    })
    
    
    return(list(
      raw = reactive(rv$raw),
      norm = reactive(rv$norm),
      final = reactive(rv$final)
    ))
  })
}

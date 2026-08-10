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
        h5("1. Name the project and author"),
        p("Enter a project name and author name. Both are carried into the generated bulk RNA-seq report."),
        h5("2. Upload counts"),
        p("Upload a raw counts file separated by comma, semicolon, or tab."),
        h5("3. Confirm detected gene ID column"),
        p("CoTRA identifies the most likely gene ID column. Check it before confirming."),
        h5("4. Choose normalization method"),
        p("Choose None, FPKM, RPKM, TPM, or TMM depending on your analysis plan."),
        h5("5. Load and confirm data"),
        p("Press Load data, review the detected gene ID column, then confirm it."),
        h5("6. Rename samples if needed"),
        p("Edit the New sample name column and click Apply sample names. Sample names must be non-empty and unique."),
        h5("7. Continue to grouping"),
        p("Use the final sample names when defining experimental groups and running downstream analysis.")
      )
    ),
    
    br(),
    
    textInput(
      ns("project_name"),
      "Project name",
      value = "",
      placeholder = "e.g. rd10 retina bulk RNA-seq"
    ),
    helpText("This project name will be included automatically in the generated report."),
    
    textInput(
      ns("author_name"),
      "Author",
      value = "CoTRA user",
      placeholder = "e.g. Umair Seemab"
    ),
    helpText("The author name will be included automatically in the generated report."),
    
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
    
    uiOutput(ns("sample_rename_ui")),
    
    strong("Current sample names"),
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
      out <- suppressWarnings(apply(mat, 2, function(x) as.numeric(x)))
      if (is.null(dim(out))) {
        out <- matrix(out, ncol = 1)
        colnames(out) <- colnames(mat)[1]
      }
      out
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
    
    apply_sample_names <- function(new_names){
      new_names <- trimws(as.character(new_names))
      
      if (any(!nzchar(new_names))) {
        showNotification("Sample names cannot be empty.", type = "error")
        return(FALSE)
      }
      
      if (anyDuplicated(new_names) > 0) {
        duplicates <- unique(new_names[duplicated(new_names)])
        showNotification(
          paste0("Sample names must be unique. Duplicated: ", paste(duplicates, collapse = ", ")),
          type = "error",
          duration = 8
        )
        return(FALSE)
      }
      
      for (nm in c("raw", "norm", "final")) {
        mat <- rv[[nm]]
        if (!is.null(mat)) {
          if (ncol(mat) != length(new_names)) {
            showNotification("Could not rename samples because the matrix dimensions changed.", type = "error")
            return(FALSE)
          }
          colnames(mat) <- new_names
          rv[[nm]] <- mat
        }
      }
      
      rv$sample_names$New <- new_names
      TRUE
    }
    
    project_name_value <- reactive({
      x <- input$project_name
      if (is.null(x)) x <- ""
      x <- trimws(x)
      if (nzchar(x)) x else "Untitled bulk RNA-seq project"
    })
    
    author_name_value <- reactive({
      x <- input$author_name
      if (is.null(x)) x <- ""
      x <- trimws(x)
      if (nzchar(x)) x else "CoTRA user"
    })
    
    rv <- reactiveValues(
      raw = NULL,
      norm = NULL,
      final = NULL,
      df = NULL,
      detected_col = NULL,
      sample_names = NULL
    )
    
    observeEvent(input$matrix_file, {
      req(input$matrix_file$name)
      current_name <- input$project_name
      if (is.null(current_name)) current_name <- ""
      
      if (!nzchar(trimws(current_name))) {
        suggested_name <- tools::file_path_sans_ext(basename(input$matrix_file$name))
        updateTextInput(session, "project_name", value = suggested_name)
      }
    }, ignoreInit = TRUE)
    
    observeEvent(input$load_btn, {
      req(input$matrix_file)
      
      rv$raw <- NULL
      rv$norm <- NULL
      rv$final <- NULL
      rv$sample_names <- NULL
      
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
      rv$sample_names <- data.frame(
        Original = colnames(norm_m),
        New = colnames(norm_m),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      
      output$id_detect_ui <- renderUI(NULL)
    })
    
    output$summary <- renderText({
      req(rv$final)
      paste0(
        "Project: ", project_name_value(),
        "\nAuthor: ", author_name_value(),
        "\nGenes: ", nrow(rv$final),
        "\nSamples: ", ncol(rv$final),
        "\nNormalization: ", input$norm_method
      )
    })
    
    output$sample_rename_ui <- renderUI({
      req(rv$sample_names)
      
      rename_rows <- lapply(seq_len(nrow(rv$sample_names)), function(i) {
        fluidRow(
          column(
            6,
            div(
              style = "padding-top:7px; word-break:break-word;",
              as.character(rv$sample_names$Original[i])
            )
          ),
          column(
            6,
            textInput(
              ns(paste0("sample_name_", i)),
              label = NULL,
              value = as.character(rv$sample_names$New[i]),
              width = "100%"
            )
          )
        )
      })
      
      tagList(
        hr(),
        h4("Rename samples"),
        helpText("Edit the new sample names below. Apply the names before defining experimental groups."),
        fluidRow(
          column(6, strong("Original sample name")),
          column(6, strong("New sample name"))
        ),
        div(
          style = "max-height:420px; overflow-y:auto; overflow-x:hidden; padding-right:8px; margin-top:8px;",
          rename_rows
        ),
        br(),
        actionButton(ns("apply_sample_names"), "Apply sample names", class = "btn btn-success"),
        actionButton(ns("reset_sample_names"), "Reset to original names", class = "btn btn-default")
      )
    })
    
    observeEvent(input$apply_sample_names, {
      req(rv$sample_names)
      
      new_names <- vapply(seq_len(nrow(rv$sample_names)), function(i) {
        value <- input[[paste0("sample_name_", i)]]
        if (is.null(value)) as.character(rv$sample_names$New[i]) else as.character(value)
      }, character(1))
      
      if (apply_sample_names(new_names)) {
        showNotification("Sample names updated. Downstream modules now use the new names.", type = "message")
      }
    })
    
    observeEvent(input$reset_sample_names, {
      req(rv$sample_names)
      original_names <- as.character(rv$sample_names$Original)
      rv$sample_names$New <- original_names
      
      for (i in seq_along(original_names)) {
        updateTextInput(session, paste0("sample_name_", i), value = original_names[i])
      }
      
      if (apply_sample_names(original_names)) {
        showNotification("Original sample names restored.", type = "message")
      }
    })
    
    output$sample_preview <- renderTable({
      req(rv$final)
      head(data.frame(Sample = colnames(rv$final), check.names = FALSE), 10)
    })
    
    output$gene_preview <- renderTable({
      req(rv$final)
      head(data.frame(GeneID = rownames(rv$final), check.names = FALSE), 10)
    })
    
    output$group_suggest <- renderText({
      req(rv$final)
      guess_groups(colnames(rv$final))
    })
    
    output$download_ui <- renderUI({
      req(rv$norm)
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
      filename = function(){ cotra_file_name("Bulk_NormalizedCountMatrix", "csv") },
      content = function(file){
        req(rv$norm)
        write.csv(rv$norm, file, row.names = TRUE)
      }
    )
    
    return(list(
      raw = reactive(rv$raw),
      norm = reactive(rv$norm),
      final = reactive(rv$final),
      project_name = project_name_value,
      author_name = author_name_value,
      sample_names = reactive(rv$sample_names)
    ))
  })
}

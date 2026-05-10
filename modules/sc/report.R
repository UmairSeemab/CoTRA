# ==========================================================
# scRNA report generator module
# ==========================================================

mod_sc_report_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    bs4Card(
      title = "scRNA-seq Report Generator",
      status = "primary",
      width = 12,
      solidHeader = TRUE,
      collapsible = TRUE,
      fluidRow(
        column(
          width = 4,
          radioButtons(
            ns("report_format"),
            "Choose format",
            choices = c("HTML" = "html", "PDF" = "pdf")
          )
        ),
        column(
          width = 4,
          textInput(
            ns("report_title"),
            "Report title",
            value = "CoTRA scRNA-seq Analysis Report"
          )
        ),
        column(
          width = 4,
          actionButton(
            ns("generate_report"),
            "Generate Report",
            class = "btn btn-success btn-block"
          )
        )
      ),
      hr(),
      downloadButton(
        ns("download_report"),
        "Download report",
        class = "btn btn-primary"
      )
    )
  )
}


mod_sc_report_server <- function(id, seurat_in, markers_in) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    report_file <- reactiveVal(NULL)
    
    observeEvent(input$generate_report, {
      req(seurat_in())
      
      showNotification("Building report", type = "message")
      
      tmp <- tempfile(fileext = ifelse(input$report_format == "pdf", ".pdf", ".html"))
      report_file(tmp)
      
      rmarkdown::render(
        input = "modules/sc/report_template.Rmd",
        output_file = tmp,
        params = list(
          title = input$report_title,
          seu = seurat_in(),
          markers = markers_in()
        ),
        envir = new.env(parent = globalenv()),
        quiet = TRUE
      )
      
      showNotification("Report ready", type = "message")
    })
    
    output$download_report <- downloadHandler(
      filename = function() {
        if (input$report_format == "pdf") "scRNA_report.pdf" else "scRNA_report.html"
      },
      content = function(file) {
        file.copy(report_file(), file)
      }
    )
  })
}

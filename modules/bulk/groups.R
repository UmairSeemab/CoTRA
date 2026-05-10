mod_bulk_groups_ui <- function(id){
  ns <- NS(id)
  
  tagList(
    shinyjs::useShinyjs(),
    
    actionButton(ns("help_toggle"), icon("question-circle"), class = "btn btn-info btn-sm"),
    
    shinyjs::hidden(
      div(
        id = ns("help_box"),
        style = "border:1px solid #d9d9d9; padding:12px; margin-top:10px; border-radius:6px; background:#f7f7f7;",
        
        h4("Step-by-Step Guide"),
        h5("1. Review uploaded samples"),
        p("The module reads sample names from the uploaded expression matrix."),
        
        h5("2. Select samples"),
        p("Choose one or more samples that belong to the same biological or experimental condition."),
        
        h5("3. Enter group name"),
        p("Type a clear group name such as Control, Treated, WT, or KO."),
        
        h5("4. Add group"),
        p("Click Add group to save the selected samples under the chosen group name."),
        
        h5("5. Repeat for all groups"),
        p("Repeat the same process until all samples are assigned to their respective groups."),
        
        h5("6. Exclude ungrouped samples"),
        p("By default, samples not assigned to any group are excluded from downstream analysis, including QC visualization."),
        
        h5("7. Review excluded samples"),
        p("The excluded samples table shows samples that are currently not assigned to any group and will be removed from downstream analysis when exclusion is enabled."),
        
        h5("8. Review or clear"),
        p("Check the group table below. Use Clear all if you want to reset the grouping and start again.")
      )
    ),
    
    br(),
    
    h4("Define experimental groups"),
    helpText("Select samples from the uploaded matrix and assign them to groups."),
    
    uiOutput(ns("sample_selector")),
    
    br(),
    
    fluidRow(
      column(6, textInput(ns("group_name"), "Group name", "")),
      column(6, actionButton(ns("add_group"), "Add group", class = "btn btn-primary"))
    ),
    
    br(),
    
    checkboxInput(
      ns("exclude_ungrouped"),
      "Exclude samples not assigned to any group from downstream analysis",
      value = TRUE
    ),
    
    hr(),
    
    DT::DTOutput(ns("group_table")),
    
    hr(),
    
    h5("Excluded samples"),
    uiOutput(ns("excluded_ui")),
    tableOutput(ns("excluded_table")),
    
    br(),
    
    actionButton(ns("clear_groups"), "Clear all", class = "btn btn-danger")
  )
}

mod_bulk_groups_server <- function(id, bulk_data){
  moduleServer(id, function(input, output, session){
    ns <- session$ns
    library(DT)
    library(shinyWidgets)
    library(shinyjs)
    
    rv <- reactiveValues(
      groups = data.frame(
        Group = character(),
        Samples = character(),
        stringsAsFactors = FALSE
      )
    )
    
    observeEvent(input$help_toggle, {
      shinyjs::toggle(id = "help_box", anim = TRUE)
    })
    
    output$sample_selector <- renderUI({
      mat <- bulk_data$final()
      req(mat)
      samples <- colnames(mat)
      
      pickerInput(
        ns("selected_samples"),
        "Select samples",
        choices = samples,
        multiple = TRUE,
        options = list(`actions-box` = TRUE)
      )
    })
    
    observeEvent(input$add_group, {
      req(input$group_name)
      req(input$selected_samples)
      
      new_row <- data.frame(
        Group = trimws(input$group_name),
        Samples = paste(input$selected_samples, collapse = "; "),
        stringsAsFactors = FALSE
      )
      
      rv$groups <- rbind(rv$groups, new_row)
    })
    
    observeEvent(input$clear_groups, {
      rv$groups <- rv$groups[0, , drop = FALSE]
    })
    
    output$group_table <- DT::renderDT({
      datatable(
        rv$groups,
        rownames = FALSE,
        options = list(scrollX = TRUE)
      )
    })
    
    grouped_samples <- reactive({
      if (nrow(rv$groups) == 0) {
        return(character(0))
      }
      
      unique(trimws(unlist(strsplit(rv$groups$Samples, ";"))))
    })
    
    ungrouped_samples <- reactive({
      mat <- bulk_data$final()
      req(mat)
      setdiff(colnames(mat), grouped_samples())
    })
    
    filtered_matrix <- reactive({
      mat <- bulk_data$final()
      req(mat)
      
      if (!isTRUE(input$exclude_ungrouped)) {
        return(mat)
      }
      
      keep_samples <- intersect(colnames(mat), grouped_samples())
      
      if (length(keep_samples) == 0) {
        return(mat[, 0, drop = FALSE])
      }
      
      mat[, keep_samples, drop = FALSE]
    })
    
    output$excluded_ui <- renderUI({
      mat <- bulk_data$final()
      req(mat)
      
      if (!isTRUE(input$exclude_ungrouped)) {
        return(tags$div("Exclusion is disabled. All samples are included."))
      }
      
      ex <- ungrouped_samples()
      
      if (length(ex) == 0) {
        return(tags$div("No excluded samples. All samples are assigned to groups."))
      }
      
      tags$div(
        style = "color:#b22222;",
        paste0("Total excluded samples: ", length(ex))
      )
    })
    
    output$excluded_table <- renderTable({
      mat <- bulk_data$final()
      req(mat)
      
      if (!isTRUE(input$exclude_ungrouped)) {
        return(NULL)
      }
      
      ex <- ungrouped_samples()
      
      if (length(ex) == 0) {
        return(NULL)
      }
      
      data.frame(
        Sample = ex,
        Status = "Excluded",
        stringsAsFactors = FALSE
      )
    })
    
    return(list(
      groups = reactive(rv$groups),
      grouped_samples = grouped_samples,
      ungrouped_samples = ungrouped_samples,
      exclude_ungrouped = reactive(input$exclude_ungrouped),
      final = filtered_matrix
    ))
  })
}
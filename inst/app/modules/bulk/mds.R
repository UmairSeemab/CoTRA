mod_bulk_mds_ui <- function(id){
  ns <- NS(id)
  plotly::plotlyOutput(ns("mds_plot"), height = "420px")
}

mod_bulk_mds_server <- function(id, bulk_data){
  moduleServer(id, function(input, output, session){
    
    output$mds_plot <- plotly::renderPlotly({
      
      req(bulk_data$final())
      
      mat <- bulk_data$final()
      
      d <- dist(t(mat))
      fit <- cmdscale(d, eig = TRUE, k = 2)
      
      df <- data.frame(
        x = fit$points[,1],
        y = fit$points[,2],
        sample = colnames(mat)
      )
      
      plotly::plot_ly(
        df,
        x = ~x,
        y = ~y,
        text = ~sample,
        type = "scatter",
        mode = "markers",
        marker = list(size = 10)
      ) |>
        plotly::layout(
          xaxis = list(title = "MDS 1"),
          yaxis = list(title = "MDS 2"),
          title = "Interactive MDS Plot"
        )
    })
  })
}

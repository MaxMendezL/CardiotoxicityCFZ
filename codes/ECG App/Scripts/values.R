HeartRateUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    selectInput(
      ns("Input"),
      "Calculate Heart Rate",
      choices = c("Control", "BTZ", "CFZ", "ATRA", "CFZATRA")
    ),
    DT::dataTableOutput(ns("values"))
  )
}


HeartRate <- function(input, output, session) {
  
  output$values <- DT::renderDataTable({
    req(input$Input)
    
    res <- analyze_ecg_signal(input$Input)
    dat <- res$table
    
    DT::datatable(
      dat,
      editable = FALSE,
      rownames = FALSE,
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        autoWidth = TRUE,
        dom = "tip"
      )
    )
  })
}


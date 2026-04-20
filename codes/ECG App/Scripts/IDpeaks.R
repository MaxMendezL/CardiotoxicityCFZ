PlotUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    selectInput(
      ns("Input"),
      "Compare to",
      choices = c("Control", "BTZ", "CFZ", "ATRA", "CFZATRA")
    ),
    plotOutput(ns("plot2"))
  )
}


customPlot2 <- function(input, output, session) {
  
  output$plot2 <- renderPlot({
    
    req(input$Input)
    
    dat <- ecg_data[[input$Input]]$raw
    
    if (!all(c("V1", "V2") %in% names(dat))) {
      stop(sprintf("Raw dataset '%s' must contain columns V1 and V2.", input$Input), call. = FALSE)
    }
    
    plot(
      dat$V1,
      dat$V2,
      type = "l",
      col = "black",
      xlab = "Time (Sec)",
      xlim = c(29, 30),
      ylab = "Voltage",
      ylim = c(-5, 5),
      main = "Raw Input",
      bty = "n"
    )
  })
  
  return()
}

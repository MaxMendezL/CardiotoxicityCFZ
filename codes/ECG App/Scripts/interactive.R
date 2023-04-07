
customPlotUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    selectInput(ns("Input"), "Data Input", choices = c("Control", "BTZ", "CFZ", "ATRA", "CFZATRA")),
    plotOutput(ns("plot"))
  )
  
}


customPlot <- function(input, output, session) {
  
  output$plot <- renderPlot ({ 
    
if (input$Input == "Control") {
  plot(Control$V1, Control$V2, type="l", col="black", xlab = "Time (Sec)",xlim=c(29,30), 
       ylab="Voltage", ylim=c(-5,5),main="Raw Input",adj=0, bty="n")
}else
  if (input$Input == "BTZ") {
    plot(BTZ$V1, BTZ$V2, type="l", col="black", xlab = "Time (Sec)",xlim=c(29,30), 
         ylab="Voltage", ylim=c(-5,5),main="Raw Input",adj=0, bty="n")
  }else
    if (input$Input == "CFZ") {
      plot(CFZ$V1, CFZ$V2,  type="l", col="black", xlab = "Time (Sec)",xlim=c(29,30), 
           ylab="Voltage", ylim=c(-5,5),main="Raw Input",adj=0, bty="n")
    }else
      if (input$Input == "ATRA") {
        plot(ATRA$V1, ATRA$V2, type="l", col="black", xlab = "Time (Sec)",xlim=c(29,30), 
             ylab="Voltage", ylim=c(-5,5),main="Raw Input",adj=0, bty="n")
      }else {
        plot(CFZATRA$V1, CFZATRA$V2, type="l", col="black", xlab = "Time (Sec)",xlim=c(29,30), 
             ylab="Voltage", ylim=c(-5,5),main="Raw Input",adj=0, bty="n")
      }
  })
  return()
}

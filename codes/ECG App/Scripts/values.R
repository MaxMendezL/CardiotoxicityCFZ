HeartRateUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    selectInput(ns("Input"), "Calculate Heart Rate", choices = c("Control", "BTZ", "CFZ", "ATRA", "CFZATRA")),
    dataTableOutput(ns("values"))
  )
}

HeartRate <- function(input, output, session) {
  
output$values <- renderDataTable(
  { 
  if (input$Input  == "Control")       {datatable(ECGvals_CTRL, editable = FALSE)   
  }else if (input$Input  == "BTZ")     {datatable(ECGvals_BTZ, editable = FALSE)
  }else if (input$Input  == "CFZ")     {datatable(ECGvals_CFZ, editable = FALSE)
  }else if (input$Input  == "ATRA")    {datatable(ECGvals_ATRA, editable = FALSE)
  }else{ 
    datatable(ECGvals_CFZATRA, editable = FALSE)
   
  }
  }

)

}



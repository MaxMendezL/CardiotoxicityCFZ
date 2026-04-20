library(shiny)
library(ggplot2)
library(dplyr)
library(pracma)
library(zoo)
library(DT)
library(shinybusy)
library(shinythemes)
library(shinyFeedback)

# Load data/environment first
source("global_env.R")

# Load modules
source("interactive.R")
source("IDpeaks.R")
source("values.R")
# source("peaks1.R")
source("Filtered.R")
source("smooth.R")
source("analysis_ecg.R")


ui <- navbarPage(
  title = "Shiny ECG",
  theme = shinythemes::shinytheme("sandstone"),
  
  tabPanel(
    "Home",
    h3("Here you can visualize the murine electrocardiograms recorded under anesthesia."),
    h4(
      "The ECGs were recorded under anesthesia. The program takes CSV files as input, ",
      "cleans background noise, identifies P and R waves, filters them, and calculates ",
      "heart rate (beats/min) and ECG intervals. All values are calculated automatically."
    )
  ),
  
  tabPanel(
    "Visualization",
    shinyFeedback::useShinyFeedback(),
    shinybusy::add_busy_spinner(spin = "fading-circle"),
    
    fluidRow(
      column(6, customPlotUI("plot")),
      column(6, PlotUI("plot2"))
    ),
    
    br(),
    br(),
    
    fluidRow(
      # column(6, PeaksUI("peaks1")),
      column(6, wavesUI("plot3")),
      column(6, HeartRateUI("values"))
    )
  )
)


server <- function(input, output, session) {
  
  callModule(customPlot, "plot")
  callModule(customPlot2, "plot2")
  # callModule(Peaks1, "peaks1")
  callModule(customPlot3, "plot3")
  callModule(HeartRate, "values")
}


shinyApp(ui = ui, server = server)   

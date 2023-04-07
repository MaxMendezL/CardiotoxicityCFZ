library(shiny)
library(ggplot2)
library(dplyr)
library(pracma)
library(zoo)
library(DT)
library(shinybusy)
library(shinythemes)
library(ggplot2)
library(shinyFeedback)

source("interactive.R")
source("IDpeaks.R")
source("values.R")
#source("peaks1.R")
source("Filtered.R")
source("global_env.R")
source("smooth.R")


ui <- navbarPage(
  "Shiny ECG",
  theme = shinytheme("sandstone"),
  tabPanel("Home",
  h3("Here you can visualize the murine electrocardiograms recorded under anesthesia as published in (our ref...)"), 
  #h4("Select in Data Input either the Control or an example of Experimental Group for visualization"),
  h4("The ECGs were recorded under anesthesia.
      The program takes csv files as input, cleans the background noise, identify P and R waves, filter them and calculates the Heart Rate (Beats/min) and ECG intervals. 
      All the values are automatically calculated")),
  tabPanel(
    "Visualization",
    fluidRow(
      useShinyFeedback(),
      add_busy_spinner(spin="fading-circle"),
      column(6, customPlotUI("plot")), 
      column(6, PlotUI("plot2"))),
    br(), 
    br(), 
    br(), 
    br(), 
    br(), 
    
    fluidRow(
      #column(6, PeaksUI("peaks1")),
      column(6, wavesUI("plot3")), 
      column(6, HeartRateUI("values")))
  ))

  

server <- function(input, output, session) {
  
    callModule(customPlot, "plot")
    callModule(customPlot2, "plot2")
    #callModule(Peaks1, "peaks1")
    callModule(customPlot3, "plot3")
    callModule(HeartRate, "values")
    
  }


shinyApp(ui=ui, server=server)   
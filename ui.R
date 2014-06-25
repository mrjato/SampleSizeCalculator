
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
# 
# http://www.rstudio.com/shiny/
#

library(shiny)

#dataTable <- function(inputId, conditions = 2) {
#  tagList(
#    singleton(tags$head(tags$script(src = "js/dataTable.js")))
#  )
#}

shinyUI(fluidPage(
  # Application title
  titlePanel("MALDI Sample Size Calculation (Presence/Absence)"),
  
  # Sidebar with a slider input for number of observations
  sidebarLayout(position = "left", fluid=TRUE,
    sidebarPanel(
      tags$h3("Parameters"),
      fluidRow(
        column(4, numericInput("numConditions", "Conditions", value="2", min=2, max=100, step=1)),
        column(7, offset=1,numericInput("numOfPeaks", "Peaks in Experiment", value=2000, min=1, max=100000, step=1))
      ),
      sliderInput("alpha", "Alpha", value=0.05, min=0.01, max=1, step=0.01),
      sliderInput("power", "Power", value=0.8, min=0.5, max=0.95, step=0.01),
      
      tags$hr(),
      tags$h3("Presence Patterns"),
      selectInput("pattern", "Pattern", 
        choices=list("1 vs 1" = "halfvshalf")
      ),
      sliderInput("presence", "Presence means", value=0.8, min=0, max=1, step=0.01, format="0%"),
      sliderInput("absence", "Absence means", value=0.2, min=0, max=1, step=0.01, format="0%"),
      conditionalPanel(condition = "input.pattern == 'onevsone'",
        sliderInput("noise", "Noise Means", value=0.5, min=0, max=1, step=0.01, format="0%")
      ),
      conditionalPanel(condition = "input.pattern == 'split'",
        uiOutput("splitLocation")
      ),
      actionButton("applyPattern", "Apply Pattern")
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      fluidRow(
        column(6, 
            tags$h3("Peak Presence"),
            uiOutput("presenceByCondition")
        ),
        column(5, offset=1,
          tags$h3("Sample Size"),
          textOutput(outputId="sampleSize"),
          textOutput(outputId="effectSize"),
          tags$hr(),
          tags$h3("Method Description"),
          htmlOutput(outputId="methodDescription")
        )
      )
    )
  )
))

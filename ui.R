# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
# 
# http://www.rstudio.com/shiny/
#

library(shiny)

shinyUI(fluidPage(
  # Application title
  titlePanel("MALDI Sample Size Calculation (Presence/Absence)"),
  
  # Sidebar with a slider input for number of observations
  tabsetPanel(
    tabPanel("Biomarker Discovery", 
      sidebarLayout(position = "left", fluid=TRUE,
        sidebarPanel(
          tags$h3("Parameters"),
          fluidRow(
            column(4, numericInput("numConditions", "Conditions", value=4, min=2, max=100, step=1)),
            column(7, offset=1,numericInput("numOfPeaks", "Peaks in Experiment", value=2000, min=1, max=100000, step=1))
          ),
          sliderInput("alpha", "Alpha", value=0.05, min=0.01, max=1, step=0.01),
          tags$p(style="font-size: 0.8em; color: gray;", 
            "Probability of detecting a peak as a potential biomarker when it is not
            a real biomarker (Type I error prob.)."
          ),
          sliderInput("power", "Power", value=0.8, min=0.5, max=0.95, step=0.01),
          tags$p(style="font-size: 0.8em; color: gray;", 
            "Probability of detecting a peak as a potential biomarker when it is is
            a real biomarker (1-beta is the Type II error prob.)."
          ),
          
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
            column(7, 
              tags$h3("Presence Pattern"),
              tags$p(style="font-size: 0.9em; color: gray;", 
                "This presence pattern allows you to easily establish your effect
                size."
              ),
              uiOutput("presenceByCondition"),
              htmlOutput("effectSize"),
              tags$hr(),
              tags$h3("Sample Heatmap (Simulation)"),
              htmlOutput("heatmapSubtitle"),
              fluidRow(
                column(6, numericInput("heatmapSamples", label="Samples by Condition", min=1, max=50, value=10, step=1)),
                column(6, numericInput("heatmapPeaks", label="Peaks", min=1, max=100, value=20, step=1))
              ),
              checkboxInput("inversePattern", label="Show also inverse pattern", value=FALSE),
              plotOutput(outputId="heatmap")
            ),
            column(5,
              tags$div(class="well",
                tags$h3("Sample Size"),
                htmlOutput(outputId="sampleSize")
              ),
              tags$hr(),
              tags$h3("Method Description"),
              htmlOutput(outputId="methodDescription")
            )
          )
        )
      )
    ),
    tabPanel("Diagnostic Test Validation", 
      tabsetPanel(
        tabPanel("Two Conditions",
         sidebarLayout(position = "left", fluid=TRUE,
           sidebarPanel(
              tags$h3("Parameters"),
              sliderInput("twoSensitivity", label="Sensitivity", min=0, max=0.99, step=0.01, value=0.9, format="0%"),
              tags$p(style="font-size: 0.8em; color: gray;", 
                "The minimum acceptable sensitivity for your diagnostic test (be conservative).
                Sensitivity is the probability of a diseased sample being diagnosted as such."
              ),
              sliderInput("twoSpecificity", label="Specificity", min=0, max=0.99, step=0.01, value=0.9, format="0%"),
              tags$p(style="font-size: 0.8em; color: gray;", 
                "The minimum acceptable specificity for your diagnostic test (be conservative).
                Specificity is the probability of a healthy sample being diagnosted as such."
              ),
              sliderInput("twoError", label="Error", min=0.01, max=0.5, step=0.01, value=0.05, format="0%"),
              tags$p(style="font-size: 0.8em; color: gray;", 
                "The maximum acceptable error for your diagnostic test (be conservative).
                It determines the width of the confidence interval. Typically 5%."
              ),
              sliderInput("twoCI", label="Confidence Inteval", min=0.8, max=0.99, step=0.01, value=0.95, format="0%"),
              tags$p(style="font-size: 0.8em; color: gray;", 
                "Probability of the accuracy fall in (accuracy-error, accuracy+error) interval.
                Typically 90%, 95% or 99%."
              ),
              selectInput("twoMode", label="Mode", 
                choices=c(
                  "Sample condition is known" = "known", 
                  "Only prevalence is known" = "prevalence"
                ),
                selected="known"
              ),
              conditionalPanel(condition="input.twoMode == 'prevalence'",
                sliderInput("twoPrevalence", label="Prevalence", min=0.01, max=0.99, step=0.01, value=0.25, format="0%"),
                tags$p(style="font-size: 0.8em; color: gray;", 
                  "Probability of being diseased in the population."
                )
              )
           ),
           mainPanel(
             tags$div(class="well well-lg",
               tags$h3("Sample Size"),
               htmlOutput("twoSampleSize")
             ),
             tags$h3("Method Description"),
             htmlOutput("twoMethodDescription")
           )
         )
        ),
        tabPanel("Two or More Conditions",
          sidebarLayout(position = "left", fluid=TRUE,
            sidebarPanel(
              tags$h3("Parameters"),
              sliderInput("tomAccuracy", label="Accuracy", min=0, max=0.99, step=0.01, value=0.8, format="0%"),
              tags$p(style="font-size: 0.8em; color: gray;", 
                "The minimum acceptable accuracy for your diagnostic test (be conservative)."
              ),
              sliderInput("tomError", label="Error", min=0.01, max=0.5, step=0.01, value=0.05, format="0%"),
              tags$p(style="font-size: 0.8em; color: gray;", 
                "The maximum acceptable error for your diagnostic test (be conservative).
                It determines the width of the confidence interval. Typically 5%."
              ),
              sliderInput("tomCI", label="Confidence Inteval", min=0.8, max=0.99, step=0.01, value=0.95, format="0%"),
              tags$p(style="font-size: 0.8em; color: gray;", 
                "Probability of the accuracy fall in (accuracy-error, accuracy+error) interval.
                Typically 90%, 95% or 99%."
              )
            ),
            mainPanel(
              tags$div(class="well well-lg",
                tags$h3("Sample Size"),
                htmlOutput("tomSampleSize")
              ),
              tags$h3("Method Description"),
              htmlOutput("tomMethodDescription")
            )
          )
        )
      )
    )
  )
))

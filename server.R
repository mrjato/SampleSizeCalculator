
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
# 
# http://www.rstudio.com/shiny/
#

library(shiny)
source('samplesize_MALDI.R')

shinyServer(function(input, output, session) {
  output$splitLocation <- renderUI({
    sliderInput("split", "Split (Presence in Conditions)", 
      value=min(input$split, input$numConditions-1),
      min=1, 
      max=input$numConditions-1, 
      step=1, 
      ticks=TRUE
    );
  })
  
  observe({ 
    if (input$numConditions == 2) {
      updateCheckboxGroupInput(session, "pattern", 
        choices=list("1 vs 1" = "halfvshalf"),
        selected="halfvshalf"
      );
    } else {
      updateCheckboxGroupInput(session, "pattern",
        choices=list(
          "Half vs Half" = "halfvshalf",
          "1 vs All" = "onevsall",
          "1 vs 1 vs Noise" = "onevsone", 
          "Split" = "split"
        ),
        selected=input$pattern
      );
    }
  });
  
  sliderValue <- reactive({switch(input$pattern,
    onevsall = function(index) {
      ifelse (index == 1, input$presence, input$absence)
    },
    onevsone = function(index) {
      ifelse (index == 1, input$presence,
              ifelse (index == 2, input$absence, input$noise)
      )
    },
    halfvshalf = function(index) {
      ifelse(index %% 2 == 1, input$presence, input$absence)
    },
    split = function(index) {
      ifelse(index <= input$split, input$presence, input$absence)
    }
  )});
  
  output$presenceByCondition <- renderUI({
    input$applyPattern;
    
    lapply(1:input$numConditions, function(index) {
      sliderInput(
        inputId=paste("samples-", index, sep=""),
        label=paste("Probability of Presence in Condition", index),
        value=sliderValue()(index),
        min=0,
        max=1,
        step=0.01,
        format="0%"
      )
    })
  });
  
  prob <- reactive({
    sapply(1:input$numConditions, function(index) {
      sampleId <- paste("samples-", index, sep="");
      x <- input[[sampleId]];
      
      c(x/input$numConditions, (1-x)/input$numConditions);
    });
  });
  
  output$sampleSize <- renderText({
    result <- computeSampleSizeMALDI(
      prob(), 
      input$numOfPeaks, 
      alpha=input$alpha, 
      power=input$power
    );
    paste("Num. of samples: ", ceiling(result$N), " (Exact: ", round(result$N, 2), ")", sep="");
  });
  
  output$effectSize <- renderText({
    result <- computeSampleSizeMALDI(prob(), input$numOfPeaks, alpha=input$alpha, power=input$power);
    
    cohenchies <- function(label) { cohen.ES(test="chisq", label)$effect.size };
    label <- ifelse (result$w<cohenchies("small"), "small", ifelse(result$w<cohenchies("medium"), "medium", "big"))
    
    paste("Effect size: ", round(result$w, 4), " (", label, ")", sep="");
  });
  
  output$methodDescription <- renderText({
    result <- computeSampleSizeMALDI(prob(), input$numOfPeaks, alpha=input$alpha, power=input$power);
    unadjustedAlpha <- input$alpha/input$numOfPeaks;
    label <- effectSizeLabel(result$w);
    
    paste(sep="",
      "<div style=\"text-align: justify\">",
        "<p>",
          "In MALDI-TOF-MS, we will consider the peak presence/absence in each 
          sample (boolean values) as potential biomarkers. &Chi;<sup>2</sup> 
          test ofindependence will be used for each candidate biomarker to 
          estimate its statistical significance (association between the 
          presence and condition (<strong>", input$numConditions, " conditions 
          tested</strong>)). In order to avoid false positives due to multiple 
          tests, and considering approximately <strong>", input$numOfPeaks, 
          " biomarkers to be tested</strong>, for an <strong>adjusted p-value 
          (with Bonferroni correction) of ", input$alpha, " (q-value)</strong>,
          we will do the power calculation for a minimum <strong>significance 
          threshold (alpha) of ", input$alpha, "/", input$numOfPeaks, "=", 
          unadjustedAlpha, " for every biomarker</strong> [1].",
        "</p>",
        
        "<p>",
          "In this sense, the power calculation for the association <strong>
          Chi-Square test with a power of ", round(input$power*100, 2), 
          "%, alpha=", unadjustedAlpha, ", and a (", label, ") &Chi;<sup>2</sup>
          effect size of ", round(result$w, 4), "</strong>, 
          <strong style=\"color: #1122AA;\">the minimal number of samples is ", 
          ceiling(result$N), "</strong>.",
        "</p>",
        
        "<p>",
          "[1] Witte JS, Elston RC, Cardon LR. On the relative sample size 
          required for multiple comparisons. Stat Med. 2000 Feb 15;19(3):369-72.
          PubMed PMID: 10649302.", 
        "</p>",
      "</div>"
    )
  });
})

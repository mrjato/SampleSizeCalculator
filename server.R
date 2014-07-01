# This file is part of Sample Size Calculator.
# 
# Copyright (c) 2014, Miguel Reboiro Jato and Daniel González Peña, 
# All rights reserved.
#
# Sample Size Calculator is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# Sample Size Calculator is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with Sample Size Calculator.If not, see <http://www.gnu.org/licenses/>.

library(shiny);
source('cohen.R');
source('samplesize_MALDI.R');
source('samplesize_gel2D.R');
source('gel2Dspots.R');

shinyServer(function(input, output, session) {
  output$splitLocation <- renderUI({
    sliderInput("split", "Split (Presence in Conditions)", 
      value=min(input$split, input$numConditions-1),
      min=1, 
      max=max(2,input$numConditions-1), 
      step=1, 
      ticks=TRUE
    );
  });
  
  observe({ 
    if (input$numConditions == 2) {
      updateCheckboxGroupInput(session, "pattern", 
        choices=list("1 vs 1" = "onevsone"),
        selected="onevsone"
      );
    } else {
      updateCheckboxGroupInput(session, "pattern",
        choices=list(
          "Half vs Half" = "halfvshalf",
          "1 vs All" = "onevsall",
          "1 vs 1 vs Noise" = "onevsonevsnoise", 
          "Split" = "split"
        ),
        selected=ifelse(input$pattern=="onevsone", "halfvshalf", input$pattern)
      );
    }
  });
  
  sliderValue <- reactive({switch(input$pattern,
    onevsall = function(index) {
      ifelse (index == 1, input$presence, input$absence)
    },
    onevsonevsnoise = function(index) {
      ifelse (index == 1, input$presence,
        ifelse (index == 2, input$absence, input$noise)
      )
    },
    onevsone = function(index) {
      ifelse(index %% 2 == 1, input$presence, input$absence)
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
        inputId=paste("condition-", index, sep=""),
        label=paste("Probability of Presence in Condition", index),
        value=sliderValue()(index),
        min=0,
        max=1,
        step=0.01,
        format="0%"
      )
    })
  });
  
  output$heatmap <- renderPlot({
    m <- matrix(nrow=input$heatmapPeaks, ncol=input$numConditions*input$heatmapSamples);
    rnames <- vector("character", nrow(m));
    cnames <- vector("character", ncol(m));
    ccolors <- vector("character", ncol(m));
    sampleColors <- cm.colors(input$numConditions);
    
    for (p in 1:nrow(m)) {
      rnames[p] <- paste("Peak", p);
      presence <- ifelse(input$inversePattern, ifelse(p <= input$heatmapPeaks/2, 1, 0), 1);
      absence <- 1 - presence;
      
      for (c in 1:input$numConditions) {
        conditionId <- paste("condition-", c, sep="");
        presenceProb <- input[[conditionId]];
        
        for (s in 1:input$heatmapSamples) {
          col <- (c-1)*input$heatmapSamples + s;
          m[p, col] <- ifelse(runif(1) <= presenceProb, presence, absence);
          cnames[col] <- paste("Sample ", c, ".", s, sep="");
          ccolors[col] <- sampleColors[c];
        }
      }
    }
    
    rownames(m) <- rnames;
    colnames(m) <- cnames;
    
    heatmap(m, 
      col=c("#00FF00FF", "#FF0000FF"),
      legend=FALSE,
      ColSideColors=ccolors,
      labRow=NA, labCol=NA,
      Rowv=NA, Colv=NA, 
      xlab="Samples (by condition)", ylab="Peaks",
      margin=c(2, 2),
      hclustfun=NA, distfun=NA, reorderfun=NA
    );
  });
  
  prob <- reactive({
    sapply(1:input$numConditions, function(index) {
      sampleId <- paste("condition-", index, sep="");
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
    paste("<strong style=\"color: #11AA22;\">", ceiling(result$N), " total samples</strong>", sep="");
  });
  
  output$effectSize <- renderText({
    result <- computeSampleSizeMALDI(prob(), input$numOfPeaks, alpha=input$alpha, power=input$power);
    
    cohenchies <- function(label) { cohen.ES(test="chisq", label)$effect.size };
    label <- ifelse (result$w<cohenchies("small"), "small", ifelse(result$w<cohenchies("medium"), "medium", "big"))
    
    paste("<strong style=\"color: #1122AA;\">Effect size: ", round(result$w, 4), " (", label, ")</strong>", sep="");
  });
  
  output$heatmapSubtitle <- renderText({
    result <- computeSampleSizeMALDI(prob(), input$numOfPeaks, alpha=input$alpha, power=input$power);
    
    cohenchies <- function(label) { cohen.ES(test="chisq", label)$effect.size };
    label <- ifelse (result$w<cohenchies("small"), "small", ifelse(result$w<cohenchies("medium"), "medium", "big"))
    
    paste("<p style=\"color: grey;\">Peaks with ", round(result$w, 4), " (", label, ") effect size.</p>", sep="");
  });
  
  output$methodDescription <- renderText({
    result <- computeSampleSizeMALDI(prob(), input$numOfPeaks, alpha=input$alpha, power=input$power);
    unadjustedAlpha <- input$alpha/input$numOfPeaks;
    label <- effectSizeLabel.chisq(result$w);
    
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
  
  output$twoSampleSize <- renderText({
    sens <- input$twoSensitivity;
    spec <- input$twoSpecificity;
    err <- input$twoError;
    
    zstar = qnorm(1 - (1 - input$twoCI) / 2);
    sensSampleSize <- zstar^2 * sens * (1-sens) / err^2;
    specSampleSize <- zstar^2 * spec * (1-spec) / err^2;
    
    if (input$twoMode == "prevalence") {
      prev <- input$twoPrevalence;
      
      if (prev < 0.5) {
        specSampleSize <- max(specSampleSize, (sensSampleSize - prev*sensSampleSize)/prev);
      } else {
        sensSampleSize <- max(sensSampleSize, (specSampleSize - (1-prev)*specSampleSize)/(1-prev));
      }
    }
    
    sensSampleSize <- ceiling(sensSampleSize);
    specSampleSize <- ceiling(specSampleSize);
    sampleSize <- sensSampleSize + specSampleSize;
    
    paste(sep="",
      "<strong style=\"color: #11AA22;\">", 
      sampleSize, 
      " total samples",
      ifelse(input$twoMode == "known", paste(" (", sensSampleSize, " diseased, ", specSampleSize, " healthy)", sep=""), ""),
      "</strong>"
    );
  });
  
  output$twoMethodDescription <- renderText({
    sens <- input$twoSensitivity*100;
    spec <- input$twoSpecificity*100;
    ci <- input$twoCI*100;
    err <- input$twoError*100;
    ciWidth <- err*2;
    
    paste(sep="",
      "<div style=\"text-align: justify;\">",
        "<p>",
          "We estimate the sample size for the minimum acceptable sensitivity (", 
          sens, "%) and specificity (", spec, "%), and a width of the ", ci, 
          "% confidence interval of ", ciWidth, "%. ",
          "Sensitivity CI: (", max(0, sens - err), "%, ", min(100, sens + err), "%). ",
          "Specificity CI: (", max(0, spec - err), "%, ", min(100, spec + err), "%). ",
          ifelse(input$twoMode == "prevalence",
            paste(sep="",
              "Since only prevalence is known, we have taken into account a prevalence of ",
              input$twoPrevalence * 100,
              "% in order to ensure that the needed number of samples for both 
              sensitivity and specificity will be obtained."
            ),
            ""
          ),
        "</p>",
      "</div>"
    )
  });
  
  output$tomSampleSize <- renderText({
    acc <- input$tomAccuracy;
    err <- input$tomError;
    
    zstar = qnorm(1 - (1 - input$tomCI) / 2);
    samplesize <- zstar^2 * acc * (1-acc) / err^2;
    
    paste("<strong style=\"color: #11AA22;\">", ceiling(samplesize), " total samples</strong>", sep="");
  });
  
  output$tomMethodDescription <- renderText({
    acc <- input$tomAccuracy*100;
    ci <- input$tomCI*100;
    err <- input$tomError*100;
    ciWidth <- err*2;
    
    paste(sep="",
      "<div style=\"text-align: justify;\">",
        "<p>",
          "We estimate the sample size for the minimum acceptable accuracy (", 
          acc, "%) and a width of the ", ci, "% confidence interval of ", 
          ciWidth, "% (", max(0, acc - err), "%, ", min(100, acc + err), "%).",
        "</p>",
      "</div>"
    )
  });
  
  output$gelSampleSize <- renderText({
    samplesize <- computeSampleSizeGel(
      input$gelCV, 
      input$gelFC, 
      input$gelPower, 
      input$gelAlpha, 
      input$gelGroups
    );
    
    paste("<strong style=\"color: #11AA22;\">", ceiling(samplesize), " samples by condition</strong>", sep="");
  });
  
  output$gelMethodDescription <- renderText({
    samplesize <- computeSampleSizeGel(
      input$gelCV, 
      input$gelFC, 
      input$gelPower, 
      input$gelAlpha, 
      input$gelGroups
    );
    testname <- testForGroups(input$gelGroups);
    effectsize <- round(computeEffectSizeGel(input$gelGroups, input$gelCV, input$gelFC), 2);
    esLabel <- effectSizeLabelForGroups(effectsize, input$gelGroups);
    power <- input$gelPower * 100;
    alpha <- input$gelAlpha;
    
    examples <- gelExamples(samplesize, input$gelPower, input$gelAlpha, input$gelGroups);
    
    examplesText <- paste(sep="", collapse="",
      "<ul>",
      apply(examples, 1, function(x) {
        paste(sep="",
          "<li>A fold-change in the mean amount of protein of at least ",
          x[1], " (", (x[1] - 1)*100, "%), unless the coefficient of variation 
          in each group is not greater than ", x[2], ".</li>"
        );
      }),
      "</ul>"
    );
    
    paste(sep="",
      "<div style=\"text-align: justify;\">",
        "<p>",
          "With a sample size of ", samplesize, " samples by condition a ", 
          testname, " for detecting statistically significant amout-of-protein 
          differences will be able to detect effect sizes of, at least, ", 
          effectsize, " (", esLabel, "), at ", power, "% of sensitivity, i.e. 
          power, and null hypothesis rejection threshold of ", alpha, ", which
          implies, for example",
          examplesText,
        "</p>",
      "</div>"
    );
  });
  
  output$gelSampleExplanation <- renderText({
    isArea <- input$gelMeasurementType == "area";
    
    paste(sep="",
      "<div style=\"text-align: justify;\">",
        "<p>",
            "The null hypothesis tested is that both conditions have the same ",
            ifelse(isArea, "area", "density"), 
            ". Therefore, the alternate hypothesis is that one of the conditions
            has a greater or lower ",
            ifelse(isArea, "area", "density"), 
            ". In this simulation we present a potential biomarker where the 
            Condition A is ",
            ifelse(isArea, "bigger", "denser"), 
            " than Condition B.<br/>",
            ifelse(isArea,
              "Stronger differences will present a minimum area in Condition A
              bigger than the maximum area in Condition B (red circles).",
              "Stronger differences will present this pattern in the same
              direction in all of the three cases (best, medium, and worst)."
            ),
        "</p>",
      "</div>"
    )
  });
  
  output$gelSampleSubtitle <- renderText({
    effectsize <- round(computeEffectSizeGel(input$gelGroups, input$gelCV, input$gelFC), 2);
    esLabel <- effectSizeLabelForGroups(effectsize, input$gelGroups);
    
    as.character(
      tags$p(style="color: grey;",
        paste(sep="", "Effect size is ", effectsize, " (", esLabel, ").")
      )
    );
  });
  
  output$gelSamplePlot <- renderPlot({
    if (input$gelMeasurementType == "area") {
      drawSpots(input$gelCV, input$gelFC);
    } else {
      drawSpotsDensity(input$gelCV, input$gelFC, baseDensity=input$gelMediumDensity);
    }
  }, height=960, width=960);
})

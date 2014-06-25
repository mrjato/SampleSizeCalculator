#Compute sample size for MALDI
library(pwr);

computeSampleSizeMALDI <-function(prob, numberofpeaks, alpha=0.05, power=0.8){
  pwr.chisq.test(
    w=ES.w2(prob),
    df=(nrow(prob)-1)*(ncol(prob)-1),
    N=NULL, 
    power=power, 
    sig.level = alpha/numberofpeaks
  );
}

cohenchies <- function(label) { cohen.ES(test="chisq", label)$effect.size };
effectSizeLabel <- function(w) { 
  ifelse (w<cohenchies("small"), "small", ifelse(w<cohenchies("medium"), "medium", "big"))
}
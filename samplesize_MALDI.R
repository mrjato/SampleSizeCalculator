#Compute sample size for MALDI
library(pwr);
source("cohen.R");

computeSampleSizeMALDI <-function(prob, numberofpeaks, alpha=0.05, power=0.8){
  pwr.chisq.test(
    w=ES.w2(prob),
    df=(nrow(prob)-1)*(ncol(prob)-1),
    N=NULL, 
    power=power, 
    sig.level = alpha/numberofpeaks
  );
}
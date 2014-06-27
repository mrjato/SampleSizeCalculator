library(pwr);

cohenD <- function(cv, fc) {
  2*(fc-1)/(cv*fc+cv)
}

cohenF <- function(cv, fc) {
  d <- cohenD(cv, fc);
  d/sqrt(d^2 + 4);
}

cohenChiES <- function(label) { 
  cohen.ES(test="chisq", label)$effect.size 
};
cohenAnovaES <- function(label) { 
  cohen.ES(test="anov", label)$effect.size 
};
cohenTtestES <- function(label) { 
  cohen.ES(test="r", label)$effect.size 
};

effectSizeLabel <- function(w, cohen) { 
  ifelse (w<cohen("small"),
    "small", 
    ifelse(w<cohen("medium"), "medium", "big")
  );
}

effectSizeLabel.chisq = function(w) {
  effectSizeLabel(w, cohenChiES);
}

effectSizeLabel.anov = function(w) {
  effectSizeLabel(w, cohenAnovaES);
}

effectSizeLabel.t = function(w) {
  effectSizeLabel(w, cohenTtestES);
}
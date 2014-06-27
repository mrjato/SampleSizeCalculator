source("cohen.R");

testForGroups <- function(groups) {
  if (groups < 2) {
    stop("groups must be >= 2");
  } else if (groups == 2) {
    "t-test";
  } else {
    "ANOVA";
  }
}

effectSizeLabelForGroups <- function(w, groups) {
  if (groups < 2) {
    stop("groups must be >= 2");
  } else if (groups == 2) {
    effectSizeLabel.t(w);
  } else {
    effectSizeLabel.chisq(w);
  }
}

computeEffectSizeGel <- function(groups, cv, fc) {
  if (groups < 2) {
    stop("groups must be >= 2");
  } else if (groups == 2) {
    cohenD(cv, fc);
  } else {
    cohenF(cv, fc);
  }
}

computeSampleSizeGel <- function(cv, fc, power, alpha, groups=2) {
  if (groups < 2) {
    stop("groups must be >= 2");
  } else if (groups == 2) {
    tryCatch({
      ceiling(pwr.t.test(
        d=cohenD(cv, fc), 
        power=power, 
        sig.level=alpha, 
        alternative="two.sided", 
        type="two.sample")$n
      );
    }, error = function(e) {
      2
    })  
  } else {
    ceiling(pwr.anova.test(
      f=cohenF(cv, fc), 
      k=groups, 
      power=power, 
      sig.level=alpha)$n
    );
  }
}

gelExamples <- function(samples, power, alpha, groups) {
  cv <- c(seq(0.1,0.6,0.1), 0.8, 1, 1.25, 1.5, 2)
  fc <- c(seq(1.1,1.6,0.1), 1.8, 2, 2.25, 2.5, 3, 4)
  
  m <- matrix(-1, nrow=length(fc), ncol=length(cv))
  rownames(m) <- fc
  colnames(m) <- cv
  
  for (f in fc) {
    fi <- as.character(f);
    for (c in cv) {
      ci <- as.character(c);
      m[fi, ci] <- computeSampleSizeGel(c, f, power, alpha, groups);
    }
  }
  
  fcs <- c();
  cvs <- c();
  
  for (f in fc) {
    fi <- as.character(f);
    bestCV <- NA
    
    for (c in cv) {
      ci <- as.character(c);
      
      if (m[fi, ci] <= samples) {
        bestCV <- c;
      }
    }
    
    if (!is.na(bestCV)) {
      fcs <- c(fcs, f);
      cvs <- c(cvs, c);
    }
  }
  
  data.frame(fc=fcs, cv=cvs);
}
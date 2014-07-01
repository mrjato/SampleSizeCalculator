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
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
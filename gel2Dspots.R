drawSpots <- function (cv, fc) {
  library(plotrix);
  
  plot(x=NULL, y=NULL,
    xlab="", ylab="", 
    xlim=c(-0.5, 0.5), ylim=c(-0.5, 0.5),
    type="n", axes=F
  );
  
  baseradius = 0.1;
  fcradius = radiusDifference(baseradius, 1 + fc);
  
  #smartzoom
  maxradius = radiusDifference(fcradius, 1 + cv);
  zoom = ifelse(maxradius > 0.25, 0.25/maxradius, 1);
  
  #applyzoom
  baseradius = baseradius * zoom;
  fcradius = fcradius * zoom;
  
  #draw fold-changed circle (1xfoldchange)
  drawSpot(-0.25, 0, fcradius, cv, desvMin="red", desvMax="green");
  text(-0.25, 0.3, labels="Condition A", cex=3);
  
  #draw reference circle (1xfoldchange)
  drawSpot(0.25, 0, baseradius, cv, desvMin="green", desvMax="red");
  text(0.25, 0.3, labels="Condition B", cex=3);
}

drawSpotsDensity <- function (cv, fc, baseDensity = 0.5) {
  library(plotrix);
  
  plot(x=NULL, y=NULL,
       xlab="", ylab="", 
       xlim=c(-1, 1), ylim=c(-1, 1),
       type="n", axes=F
  );
  
  fcAlpha = min(1, baseDensity * fc);
  
  drawSpot(0.25, -0.5, 0.2, dotCol=rgb(0, 0, 1, min(1, baseDensity * (1+cv))));
  drawSpot(0.25, 0, 0.2, dotCol=rgb(0, 0, 1, baseDensity));
  drawSpot(0.25, 0.5, 0.2, dotCol=rgb(0, 0, 1, max(0, baseDensity * (1-cv))));
  
  drawSpot(-0.25, 0.5, 0.2, dotCol=rgb(0, 0, 1, min(1, fcAlpha * (1+cv))));
  drawSpot(-0.25, 0, 0.2, dotCol=rgb(0, 0, 1, fcAlpha));
  drawSpot(-0.25, -0.5, 0.2, dotCol=rgb(0, 0, 1, max(0, fcAlpha * (1-cv))));
  
  text(-0.25, 0.8, labels="Condition A", cex=2);
  text(0.25, 0.8, labels="Condition B", cex=2);
  text(-0.5, 0.5, labels="Best", adj=1, cex=2);
  text(-0.5, 0, labels="Medium", adj=1, cex=2);
  text(-0.5, -0.5, labels="Worst", adj=1, cex=2);
}

radiusDifference <- function(radius, difference) {
  area <- pi * radius * radius;
  newarea <- area * difference;
  
  ifelse(newarea > 0, max(0, sqrt(newarea/pi)), 0);
}

drawSpot <- function(x, y, radius=1, desv=0, dotCol="black", desvMin="red", desvMax="green") {
  draw.circle(x, y, radius, border=NA, col=dotCol);
  
  if (desv > 0) {
    minDesv <- radiusDifference(radius, 1 - desv);
    maxDesv <- radiusDifference(radius, 1 + desv);
    draw.circle(x, y, minDesv, border=desvMin, col=NA, lwd=2, lty=5);
    draw.circle(x, y, maxDesv, border=desvMax, col=NA, lwd=2, lty=5);
  }
}
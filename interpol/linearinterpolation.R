linearinterpolation <- function(x1,y1,x2,y2,A) {
  %#Function requires two points (x1,x2) and their respective evaluated values at those points (y1,y2)
  %#and point A between x1 and x2 for which you want to approximate yA.
  yA=y1+(A-x1)*((y2-y1)/(x2-x1));
}
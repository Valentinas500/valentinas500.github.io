function[yA]=linearinterpolation(x1,y1,x2,y2,A)
%Function requires two points (x1,x2) and evaluated values at those points (y1,y2), and point A between x1 and x2 for which you want to approximate y(A).
yA=y1+(A-x1)*((y2-y1)/(x2-x1));
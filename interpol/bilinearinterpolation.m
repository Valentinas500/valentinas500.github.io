function[fksFinal]=bilinearinterpolation(x1,x2,y1,y2,fun,X,Y)
%Function requires endpoints of x range (x1,x2), endpoints of y range (y1,y2),
%a function (fun) with two variables,
%and point (X,Y) for which you want the function to be approximated.
%Example: bilinearinterpolation(1,10,1,10,@(x,y)(x+y)^2,5,5) 

%Evaluate the given expression at endpoints
fQ11=fun(x1,y1);
fQ12=fun(x1,y2);
fQ21=fun(x2,y1);
fQ22=fun(x2,y2);

%First, we do linear interpolation in the x-direction.
fx1= fQ11*((x2-X)/(x2-x1)) + fQ21*((X-x1)/(x2-x1));
fx2= fQ12*((x2-X)/(x2-x1)) + fQ22*((X-x1)/(x2-x1));

%Then, we interpolate in the y-direction.
fksFinal= fx1*((y2-Y)/(y2-y1)) + fx2*((Y-y1)/(y2-y1));
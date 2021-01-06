def bilinearinterpolation(X1,X2,Y1,Y2,fun,X,Y):
  #Function requires endpoints of x range (X1,X2), endpoints of y range (Y1,Y2),
  #a function (fun) with two variables x1 and x2 (using lowercase in quotations), as for example "x1+x2",
  #and point (X,Y) for which you want the function to be approximated.
  #Example: bilinearinterpolation(1,10,1,10,"pow(x1+x2,2)",5,5) 
  
  # Evaluate the given expression at endpoints
  x1 = X1; x2 = Y1; fQ11=eval(fun);
  x1 = X1; x2 = Y2; fQ12=eval(fun);
  x1 = X2; x2 = Y1; fQ21=eval(fun);
  x1 = X2; x2 = Y2; fQ22=eval(fun);
  
  #First, we do linear interpolation in the x-direction.
  fx1= fQ11*((X2-X)/(X2-X1)) + fQ21*((X-X1)/(X2-X1));
  fx2= fQ12*((X2-X)/(X2-X1)) + fQ22*((X-X1)/(X2-X1));
  
  #Then, we interpolate in the y-direction.
  fksFinal= fx1*((Y2-Y)/(Y2-Y1)) + fx2*((Y-Y1)/(Y2-Y1));
  return fksFinal

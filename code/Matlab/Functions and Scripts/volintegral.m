function [vol] = volintegral(func, xdata, ydata, matorfunc)
  % function for evaluating volume under a matrix. func must be of the form f(x, y) and be defined over xdata and ydata
  % put something in the fourth argument if func is a function handle - probably doesn't work
  if nargin<4 
    Q = func ;
  else 
    Q = func(xdata, ydata) ;
  end
  vol = trapz(ydata, trapz(xdata, Q, 2), 1) ;

function [x, fx,fcnt] = findSpeak(K,R,v,region,tol, dir);
%function [x, fx,fcnt] = findSpeak(K,R,v,mf_region,tol);
% find function (S=f(x)) peak 
%  
% Input is:     K       wavenumber
%               v       array locations (starting at 0)
%               region  [xmin xmax ymin ymax]
%               tol     frac. prec. of result
%               dir     direction in which to evaluate S
%
% Output is:    x       location of peak
%               fx      value of function S at x
%               fcnt    no. of function evals.
% Method: Golden mean sectioning, Num. Rec. routine "golden"
% S. Singh 3/02
% assumes max. is bracketed initially by end points.

% $Id: findSpeak.m 936 2004-04-23 02:29:49Z matt $

gr = 0.618;

if dir==0
  ax = region(1);
  cx  = region(2);
  y = region(3);
else
  ax = region(3);
  cx  = region(4);
  y = region(1);
end  

  bx = ((1-gr)*cx+ax*gr);

x0 = ax;
x3=cx;

if (abs(cx-bx) < abs(bx-ax))
  x1 = bx;
  x2 = bx+(1-gr)*(cx-bx);
else
  x2 = bx;
  x1 = bx-(1-gr)*(bx-ax);
end

  f1= calcS(K,R,v,x1,y,dir);
f2 = calcS(K,R,v,x2,y,dir);
fcnt=2;
while (abs(x3-x0)>tol*(abs(x1)+abs(x2)))
  if (f2>f1)
    x0=x1;x1=x2;x2=gr*x1+(1-gr)*x3;
    f1=f2;f2=calcS(K,R,v,x2,y,dir);
  else
    x3=x2;x2=x1;x1=gr*x2+(1-gr)*x0;
    f2=f1;f1=calcS(K,R,v,x1,y,dir);
  end
  fcnt = fcnt+1;
end

if (f1>f2)
  x = x1;
  fx = f1;
else
  x = x2;
  fx = f2;
end


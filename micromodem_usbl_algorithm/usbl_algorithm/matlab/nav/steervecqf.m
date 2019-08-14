% Make a steering vector set and evaluate S along it
%
% Input is:     K       wavenumber
%               v       array locations (starting at 0)
%               region  [xmin xmax ymin ymax]
%               res     resolution of search
%
% Output is:    A       steering vector
%               xvec    x axis for the P matrix
%               yvec    y axis for the P matrix
% [A, xvec, yvec] = steervec(K, v, region, res)

% $Id: steervecqf.m 936 2004-04-23 02:29:49Z matt $

function [A,S] = steervecqf (k,R, v,region,res)

    ax = region(1);
    axend = region(2);
    ayend = region(4);

		    
    xind = 1;
    acnt = 1;
    while ax <= axend

      yind = 1;
      ay = region(3);
      while ay <= ayend

	a = exp(j*k*v*[ax;ay]);
	A(:,acnt) = a;
	acnt = acnt + 1;
	S(xind,yind) = a'*R*a;
	yvec(yind) = ay;
	ay = ay + res;
	yind = yind + 1;
      end;
      xvec(xind) = ax;
      ax = ax + res;
      xind = xind + 1;
    end;



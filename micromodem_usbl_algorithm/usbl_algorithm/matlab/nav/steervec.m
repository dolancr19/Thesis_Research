% Make a steering vector set.
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

function [A, xvec, yvec] = steervec (k,v,region,res)

    ax = region(1);
    axend = region(2);
    ayend = region(4);
    A = [];
		    
    xind = 1;
    acnt = 1;
    while ax <= axend

      yind = 1;
      ay = region(3);
      while ay <= ayend

	a = exp(j*k*v*[ax;ay]);
	A(:,acnt) = a;
	acnt = acnt + 1;
	yvec(yind) = ay;
	ay = ay + res;
	yind = yind + 1;
      end;
      xvec(xind) = ax;
      ax = ax + res;
      xind = xind + 1;
    end;



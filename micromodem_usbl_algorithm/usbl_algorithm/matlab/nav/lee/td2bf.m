% Function to perform beamforming using an exhaustive signal arrival space
% search for an arbitrary 2-d array.
%
% Input is:     R       autocorrelation matrix
%               v       array locations (starting at 0)
%               f       signal frequency
%               c       speed of sound
%               region  [xmin xmax ymin ymax]
%               res     resolution of search
%
% Output is:    P       the output power
%               xvec    x axis for the P matrix
%               yvec    y axis for the P matrix
% [P, xvec, yvec] = td2bf(R, v, f, c, region, res)

function [P, xvec, yvec, A] = td2bf (R,v,f,c,region,res)


    lambda = c/f;
    k = 2*pi/lambda;

    ax = region(1);
    axend = region(2);
    ayend = region(4);

		    
    xind = 1;
    acnt = 1;
    while ax < axend

      yind = 1;
      ay = region(3);
      while ay < ayend

	a = exp(j*k*v*[ax;ay]);
	A(:,acnt) = a;
	acnt = acnt + 1;
	P(xind,yind) = (abs((a)'*R*a));
	yvec(yind) = ay;
	ay = ay + res;
	yind = yind + 1;
      end;
      xvec(xind) = ax;
      ax = ax + res;
      xind = xind + 1;

    end;



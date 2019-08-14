% Function to compute quadratic form 
%
% Input is:     R       cross-correlation matrix
%               c       speed of sound
%               region  [xmin xmax ymin ymax]
%               res     resolution of search
%
% Output is:    S       
%
% [S] = qf(R, A, region, res)

function [S] = td2bf (R,A,region,res)

    ax = region(1);
    axend = region(2);
    ayend = region(4);
    S =[];
		    
    xind = 1;
    acnt = 1;
    while ax <= axend

      yind = 1;
      ay = region(3);
      while ay <= ayend

	a = A(:,acnt);
	acnt = acnt + 1;
	S(xind,yind) = a'*R*a;
	ay = ay + res;
	yind = yind + 1;
      end;
      xvec(xind) = ax;
      xind = xind + 1;
      ax = ax + res;
    end;



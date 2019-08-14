% evaluate S at a point ax, ay
%
% Input is:     K       wavenumber
%               v       array locations (starting at 0)
%               region  [xmin xmax ymin ymax]
%               res     resolution of search
%
% Output is:    A       steering vector
%               xvec    x axis for the P matrix
%               yvec    y axis for the P matrix
% function S = calcS (k,R, v,ax,ay)

% $Id: calcS.m 936 2004-04-23 02:29:49Z matt $

function S = calcS (k,R, v,ax,ay,dir)

if (dir==1)
  tmp=ax;
  ax=ay;
  ay=tmp;
end

  a = exp(j*k*v*[ax;ay]);
S = a'*R*a;


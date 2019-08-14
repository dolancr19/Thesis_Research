% Function to perform beamforming using a golden mean minimization 
% search for an arbitrary 2-d array.
%
% Input is:     R       autocorrelation matrix
%               v       array locations (starting at 0)
%               f       signal frequency
%               c       speed of sound
%               region  [xmin xmax ymin ymax]
%               res     resolution of search
%               debug   verbosity flag  (optional)
%
% Output is:    [x_unit_val y_unit_val]
% [unitvals, az_el] = bf_fast(R, v, f, c, region, res,res_loops,debug)

% $Id: bfgolden.m 6840 2010-08-18 13:50:12Z sandipa $

function [unitvals, az_el, ncalls] = bfgolden (R,v,f,c,region,res)

if nargin < 8   debug = 0;      end;

lambda = c/f;
K = 2*pi/lambda;
array_locs = v;
mf_region = region ;            % local copy of region, will be updated below.
mf_res = res ;                  % local copy of res, will be updated below.
xmin = region(1);
xmax = region(2);
ymin = region(3);
ymax = region(4);
xline = (xmin+xmax)./2 ;
yline = (ymin+ymax)./2  ;
ncalls = 0;
maxloops = 30;
maxS1 = 0;
xprev = 0;yprev=0;
for (i=1:maxloops)

  % x line search first
  mf_region = [xmin xmax yline yline];
[xline, maxS0, nevals] = findSpeak(K,R,array_locs,mf_region,mf_res,0);
  ncalls = ncalls + nevals;
  % if function max acheived, break

%  if abs(maxS0-maxS1) < 0.0005, break; end;
  if abs(maxS0)<=abs(maxS1), xline = xprev;break; end;
xprev = xline;
  % y line search
  mf_region = [xline xline ymin ymax] ;
[yline, maxS1,nevals] = findSpeak(K,R,array_locs,mf_region,mf_res,1);
  ncalls = ncalls + nevals;
 % abs(maxS1)
%  if abs(maxS0-maxS1) < 0.05, break; end;
  if abs(maxS1)<=abs(maxS0), yline = yprev; break; end;
yprev = yline;
end;    % end loop over resolution

% change signs to fix coordinate system
xunitval = xline;
yunitval = yline;
unitvals = [ xunitval yunitval] ;
tmp = xunitval*xunitval + yunitval*yunitval ;
zunitval = sqrt(1.1 - tmp);
if 0
    if tmp <=1
        zunitval = sqrt(1 - tmp);
    else
        zunitval = 0;
        %xunitval = xunitval./sqrt(tmp);
        %yunitval = yunitval./sqrt(tmp);
    end
end

if xunitval == 0.0
  xunitval = eps;
end;
if yunitval == 0.0
  yunitval == eps;
end;

%az = atan2(yunitval,xunitval) * (180/pi);


if real(zunitval) == 0.0
  zunitval = 0.0;
end;
az = atan2(xunitval,zunitval) * (180/pi);

%el = atan2(zunitval,sqrt(tmp)) * (180/pi);
el = atan2(yunitval,zunitval) * (180/pi);
az_el = [az el];
unitvals = [xunitval, yunitval, zunitval];






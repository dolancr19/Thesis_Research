% Function to perform beamforming using an line 
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

function [unitvals, az_el] = bf_fast (R,v,f,c,region,res,res_loops,debug)

if nargin < 8   debug = 0;      end;

lambda = c/f;
K = 2*pi/lambda;
array_locs = v;
mf_region = region ;            % local copy of region, will be updated below.
mf_res = res ;                  % local copy of res, will be updated below.
res_cnt = 1 ;
xmin = region(1);
xmax = region(2);
ymin = region(3);
ymax = region(4);
xline = 0.0 ;
yline = 0.0 ;
xrange = xmax - xmin ;
yrange = ymax - ymin ;
max_x_loc = 0 ;
max_y_loc = 0 ;
rangediv = 2;
resdiv = 2;
while res_cnt <= res_loops

	% x line search first
	mf_region = [xmin xmax yline yline];
	Ax = steervec(K,array_locs,mf_region,mf_res);
	Sx = qf(R,Ax,mf_region,mf_res);


	[max_x_val,max_x_loc] = max(Sx);

	% The peak of the matrix S indicates the x-y unit vector.
	% Get the location of that peak and then get from the region and
	% resolution the actual x-y unit vector values.
	xline = xmin + (max_x_loc-1) * mf_res; 

	% y line search
	mf_region = [xline xline ymin ymax] ;
	Ay = steervec(K,array_locs,mf_region,mf_res);

	Sy = qf(R,Ay,mf_region,mf_res) ;

	[max_y_val,max_y_loc] = max(Sy) ;
	yline = ymin + (max_y_loc-1) * mf_res ;
	if debug
	  fprintf('%d res %f x-yline %f  %f x-yrange %f %f\n',res_cnt, mf_res, xline, yline, xrange, yrange);
	end;
	mf_res = mf_res/resdiv ;
	xrange = xrange/resdiv;
	yrange = yrange/resdiv;
	xmin = xline - xrange/rangediv ;
	if xmin<-1 xmin = -1; end;
	xmax = xline + xrange/rangediv ;
	if xmax>1 xmax = 1; end;
	ymin = yline - yrange/rangediv ;
	if ymin<-1 ymin = -1; end;
	ymax = yline + yrange/rangediv ;
	if ymax>1 ymax = 1; end;

	res_cnt = res_cnt + 1;
end;    % end loop over resolution

% change signs to fix coordinate system
xunitval = xline;
yunitval = yline;
unitvals = [ xunitval yunitval] ;
tmp = xunitval*xunitval + yunitval*yunitval ;
if tmp < 1.0
  zunitval = sqrt(1 - tmp);
else
  zunitval = 0.0 ;
end;
if xunitval == 0.0
  xunitval = eps;
end;
if yunitval == 0.0
  yunitval == eps;
end;

az = atan(xunitval/zunitval) * (180/pi);


if real(zunitval) == 0.0
  zunitval = 0.0;
end;

el = atan(yunitval/zunitval) * (180/pi);
az_el = [az el];



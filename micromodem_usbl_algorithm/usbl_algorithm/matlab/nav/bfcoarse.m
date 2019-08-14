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

% $Id: bfcoarse.m 6840 2010-08-18 13:50:12Z sandipa $

function [search_area] = bfcoarse(R,v,f,c,region,res)

if nargin < 8   debug = 0;      end;

lambda = c/f;
K = 2*pi/lambda;
for kk=1:size(v,2)
    if length(find(v(:,kk)==0)) == size(v,2)
        break;
    end
end
array_locs = v;
array_locs(:,kk) = [];
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
sxprev = 0;
syprev = 0;

	% x line search first
	mf_region = [xmin xmax yline yline];
	[Ax,Sx] = steervecqf(K,R,array_locs,mf_region,mf_res);


	[max_x_val,max_x_loc] = max(Sx);

	% The peak of the matrix S indicates the x-y unit vector.
	% Get the location of that peak and then get from the region and
	% resolution the actual x-y unit vector values.
	xline = xmin + (max_x_loc-1) * mf_res; 

	% y line search
	mf_region = [xline xline ymin ymax] ;
	[Ay Sy] = steervecqf(K,R,array_locs,mf_region,mf_res);
	figure(1);subplot(2,1,1); plot([xmin:mf_res:xmax-mf_res], abs(Sx));
	figure(1);subplot(2,1,2); plot([ymin:mf_res:ymax-mf_res], abs(Sy));

	[max_y_val,max_y_loc] = max(Sy) ;
        yline = ymin + (max_y_loc-1) * mf_res ;
 search_area = [max([xmin, xline-2*res]),min([xmax, xline+2*res]),...
     max([ymin, yline-2*res]),min([ymax, yline+2*res])];

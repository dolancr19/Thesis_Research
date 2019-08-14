% Function to perform beamforming using an line
% search for an arbitrary 2-d array.
%
% Input is:     R       autocorrelation matrix
%               v       array locations (starting at 0), allow 3-D
%               locations
%               f       signal frequency
%               c       speed of sound
%               region  [xmin xmax ymin ymax]
%               res     resolution of search
%               debug   verbosity flag  (optional)
%
% Output is:    [x_unit_val y_unit_val]
% [unitvals, az_el] = bf_fast(R, v, f, c, region, res,res_loops,debug)

function [unitvals, az_el, res_cnt] = bf_fast (R,v,f,c,region,res,res_loops,debug)

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
while ((res_cnt <= res_loops) & (xmin < xmax))
    
    % x line search first
    mf_region = [xmin xmax yline yline];
    Ax = steervec(K,array_locs,mf_region,mf_res);
    Sx = qf(R,Ax,mf_region,mf_res);
    
    
    [max_x_val,max_x_loc] = max(Sx);
    
    % The peak of the matrix S indicates the x-y unit vector.
    % Get the location of that peak and then get from the region and
    % resolution the actual x-y unit vector values.
    xline = xmin + (max_x_loc-1) * mf_res;
    figure(1);subplot(2,1,1); plot(abs(Sx));
    title(sprintf('%.02f x-line, %.02f y-line, %f res\n',xline, yline, mf_res));
    % y line search in a unit circle
    ylim = sqrt(1-xline.^2);
    if ymin<-ylim ymin = -ylim; end;
    if ymax>ylim ymax = ylim; end;
    mf_region = [xline xline ymin,ymax] ;
    Ay = steervec(K,array_locs,mf_region,mf_res);
    
    Sy = qf(R,Ay,mf_region,mf_res) ;
    if isempty(Sy)
        break;
    end
    [max_y_val,max_y_loc] = max(Sy) ;
    yline = ymin + (max_y_loc-1) * mf_res ;
    figure(1);subplot(2,1,2); plot(abs(Sy));
    if debug
        fprintf('%d res %f x-yline %f  %f x-yrange %f %f\n',res_cnt, mf_res, xline, yline, xrange, yrange);
    end;
    mf_res = mf_res/resdiv ;
    xrange = xrange/resdiv;
    yrange = yrange/resdiv;
    xlim = sqrt(1-yline.^2);
    ylim = sqrt(1-xline.^2);
    xmin = xline - xrange/rangediv ;
    if xmin<-xlim xmin = -xlim; end;
    xmax = xline + xrange/rangediv ;
    if xmax>xlim xmax = xlim; end;
    ymin = yline - yrange/rangediv ;
    if ymin<-ylim ymin = -ylim; end;
    ymax = yline + yrange/rangediv ;
    if ymax>ylim ymax = ylim; end;
    
    res_cnt = res_cnt + 1;
end;    % end loop over resolution

% change signs to fix coordinate system
xunitval = xline;
yunitval = yline;
tmp = xunitval*xunitval + yunitval*yunitval;
zunitval = sqrt(1.00 - tmp);
if 0
    if tmp <=1
        zunitval = sqrt(1 - tmp);
    else
        zunitval = 0;
        xunitval = xunitval./sqrt(tmp);
        yunitval = yunitval./sqrt(tmp);
    end
end
%az = atan2(xunitval,yunitval) * (180/pi);

unitvals = [ xunitval yunitval zunitval];
if kk ==1 % we were actually in the x-z plane, switch the unitvals accordingly
    zunitval = -unitvals(1);
    yunitval = unitvals(3);
    xunitval = -unitvals(2);
end
if kk ==2 % we were actually in the x-z plane, switch the unitvals accordingly
    zunitval = -unitvals(2);
    yunitval = -unitvals(3);
    xunitval = unitvals(1);
end

unitvals = [ xunitval yunitval zunitval];
az = atan2(xunitval,zunitval) * (180/pi);

%el = atan2(zunitval,sqrt(tmp)) * (180/pi);
el = atan2(yunitval,zunitval) * (180/pi);
az_el = [az el];



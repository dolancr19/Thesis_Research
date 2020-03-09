%% Helper function: dataInterp

function [cGridVal,dGridVal] = dataInterp(cmat,depmat)
% If you plan to use MATLAB's interp function, you should make an gridded
% interpolant first and then put them in a structure. The reason you do
% this is purely for computational speed. Otherwise, you could do a
% weighted average interpolation pretty quickly in the helper functions,
% getDepth and getCVal, and you don't need this function at all.

% If you'd like to discuss the differences between these two methods, reach
% out to Eeshan.

cGridVal = 0;
dGridVal = 0;

% create sound speed interpolation objects to avoid redundant calculations
depth = cmat(:,1); 
ssp = cmat(:,2);
cGridVal= griddedInterpolant(depth,ssp);

% create depth interpolation object to avoid redundant calculations
range = depmat(:,1); 
bathy = depmat(:,2);
dGridVal = griddedInterpolant(range,bathy);

end
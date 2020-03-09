% Christopher Dolan

%% Set up workspace 
clear variables
clc

% The form of the function for running the ray tracing is:
% raytrace(sstep,depmat,cmat,thetas,zsources,f0,numsteps,debug)

% sstep     the step size along the ray
% depmat    the depth matrix indexed by range, with form [r1,dep1;r2,dep2...]
% cmat      the soundspeed profile, with form [z1 c1;z2 c2...]
% thetas    the set of beam angles you want to run the raytrace code at,
%           between -90 and 90 degrees
% zsources  the source depths you want for beams

load SSP_498
load depmatdeep
load under
sstep = 10;
thetas = ele_under*pi/180;
calc_new = zeros(4*length(thetas),200);
zsources = 3;
build_calc_new(sstep,depmat,cmat,thetas,zsources,calc_new)
load calc_new


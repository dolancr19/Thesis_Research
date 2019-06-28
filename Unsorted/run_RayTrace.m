% Christopher Dolan

%% Set up workspace 
close all
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
% f0        frequency
% numsteps  number of steps to run
% debug     flag that, if true, creates plots with each top/bottom interaction


% if you want a different soundspeed profile, the format is [z1,c1;z2,c2;...]
% for example, c_munk is constructed as:
% depvec=(1:5000)';
% zbar = 2.*(z-1300)./1300;
% cvec=1500.*(1+.00737.*(zbar-1+exp(-zbar)));
% cmunk = [depvec,cvec];

% Global Parameters
sstep=10; % size of steps along ray in m
f0 = 10; % frequency
load cmunkdeep; %Munk profile to 7000m
load depmatdeep; %flat bottom, contains depmat

z_glider = 2000; %glider depth in m
theta_glider = 1.3089; %acoustic incident angle at glider in rad
tau_glider = 1.3716; %acoustic ray travel time in sec

theta_error = (3*pi/180); %error model to determine analysis range, in degrees

%calculate range of thetas for analysis based on incident angle and error
%range
theta_surf_plus = acos((cmunkdeep(1,2)/cmunkdeep(1+z_glider,2))*cos(round(theta_glider,2) + theta_error));
theta_surf_minus = acos((cmunkdeep(1,2)/cmunkdeep(1+z_glider,2))*cos(round(theta_glider,2) - theta_error));

%specify range of thets
thetas = theta_surf_minus:(theta_error*2)/30:theta_surf_plus;

%source depth, set to zero unless a source lowered into water at a
%specified depth
zsources = 0;

%call the raytrace function to calculate the ray paths based on the thetas
%chosen
raytrace(sstep,depmatdeep,cmunkdeep,thetas,zsources,f0,false,z_glider)

%call function to determine best horizontal position based on incident
%angle, depth and travel time
[r_out] = get_r(thetas, z_glider, tau_glider)

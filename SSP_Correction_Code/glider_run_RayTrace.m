% Christopher Dolan
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

% Global Parameters
sstep=10; % size of steps along ray in m
f0 = 10; % frequency
load cmunkdeep; %Munk profile to 7000m
load depmatdeep; %flat bottom, contains depmat

z_glider = 4401; %glider depth in m
theta_glider = 1.3101; %acoustic incident angle at glider in rad
tau_glider = 2.9944; %acoustic ray travel time in sec

theta_error = (3*pi/180); %error model to determine analysis range, in degrees

%calculate range of thetas for analysis based on incident angle and error
%range
theta_surf_plus = acos((cmunkdeep(1,2)/cmunkdeep(1+round(z_glider),2))*cos(theta_glider + theta_error));
theta_surf_minus = acos((cmunkdeep(1,2)/cmunkdeep(1+round(z_glider),2))*cos(theta_glider - theta_error));

%specify range of thetas
thetas = theta_surf_minus:(theta_error*2)/30:theta_surf_plus;

%source depth, set to 8m based on Liquid Robotics Waveglider literature
zsources = 8;

%call the raytrace function to calculate the ray paths based on the thetas
%chosen
glider_raytrace(sstep,depmatdeep,cmunkdeep,thetas,zsources,z_glider)

%call function to determine best horizontal position based on incident
%angle, depth and travel time
[r_out] = get_r(thetas, z_glider, tau_glider)

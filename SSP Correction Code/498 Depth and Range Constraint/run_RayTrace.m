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

sstep = 10; % size of steps along ray in m
div = 1; %number of steps between each r_out
theta_error = (3*pi/180); %error model to determine analysis range, in degrees

load SSP_498; 
load depmatdeep; %flat bottom, contains depmat
load under
calc_r = zeros(length(ele_under(1,:)),4);

for aa = 1:length(ele_under(1,:))
   % for ll = 1:length(dd)
        theta_glider = (pi/180)*ele_under(1,aa); %acoustic incident angle at glider in rad
        r_glider = hr_under(1,aa); 
        
        %tau_glider = calc_new(bb+2,dd(ll)); %acoustic ray travel time in sec
        z_glider = depth_under(1,aa);
        %calculate range of thetas for analysis based on incident angle and error
        %range
        theta_surf_plus = theta_glider + theta_error;
        theta_surf_minus = theta_glider - theta_error;

        %specify range of thetas
        thetas = theta_surf_minus:(theta_error*2)/30:theta_surf_plus;

        %source depth, set to 3m based on glider construction
        zsources = 3;

        %call the raytrace function to calculate the ray paths based on the thetas
        %chosen
        raytrace(sstep,depmatdeep,cmat,thetas,zsources,z_glider)

        %call function to determine best horizontal position based on incident
        %angle, depth and travel time
        [r_out] = get_r(thetas, z_glider, r_glider);
        calc_r(aa,:) = r_out;
        save('calc_r', 'calc_r')
    %end
end

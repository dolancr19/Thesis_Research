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
div = 10; %number of steps between each r_out
theta_error = (3*pi/180); %error model to determine analysis range, in degrees

load cmunkdeep; %Munk profile to 7000m
load depmatdeep; %flat bottom, contains depmat
load calc_new
calc_r = zeros(.25*length(calc_new(:,1)),1000);

for aa = 1:.25*length(calc_new(:,1))
    %define index to find the correct row for each angle
    bb = aa + 3*(aa-1);
    %determine the last non-zero cell in the array
    cc = 1;
    while calc_new(bb,cc) > 0
        cc = cc+1;
    end
    cc = cc-1;
    dd = div:div:cc;
    if dd(end) ~= cc
        dd = [dd cc];
    else
    end
    
    %iterate the horizontal range calculation for each stored depth
    for ll = 1:length(dd)
        z_glider = calc_new(bb,dd(ll)); %glider depth in m
        theta_glider = calc_new(bb+1,dd(ll)); %acoustic incident angle at glider in rad
        tau_glider = calc_new(bb+2,dd(ll)); %acoustic ray travel time in sec
        
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
        raytrace(sstep,depmat,cmat,thetas,zsources,z_glider)

        %call function to determine best horizontal position based on incident
        %angle, depth and travel time
        [r_out] = get_r(thetas, z_glider, tau_glider);
        calc_r(aa,ll) = r_out;
        save('calc_r', 'calc_r')
    end
end

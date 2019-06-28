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

load cmunkdeep
load depmatdeep
sstep = 10;
thetas=[90,85,80,75,70,65,60,55,50,45]*pi/180;
calc_new = zeros(4*length(thetas),1000);
zsources = 3;
%% Interpolate Data
[cVal,dVal] = dataInterp(cmat,depmat);

%% Ray Tracing

% loop over all sources
for dzs=1:length(zsources)
    z0=zsources(dzs); % get source depth
    
    % loop over all thetas
    for dth=1:length(thetas) % get theta value
        %% Initialize variables needed from raytracing
        
        % get sound speed, 1st derivative of sound speed, 2nd derivative of sound speed
        [c0,cz0,czz0]=getCVal(0,z0,cVal);
        
        % for each theta value, calculate ray
        theta0 = thetas(dth);
        
        % initialize raytrace vectors
        sray=[0];%s-values for ray
        xiray=[cos(theta0)/c0];% xi values for ray
        zetaray=[sin(theta0)/c0];% zeta values for ray
        thetavec=[theta0]; % theta values for ray
        rray=[0];%range vector for ray
        zray=[z0];%depth vector for ray
        cray=[c0];%soundspeed vector for ray
        tau=[0];%tau vector for ray
        cz=[cz0]; %first derivative of ssp
        czz=[czz0]; %second derivative of ssp
        i=1; %step count
        
        
        %% Iterate over numsteps to calculate out the ray
        
        while zray<6000
            
            r=rray(i);%r from previous step
            z=zray(i);%z from previous step 
                      
            r = r + sstep*cray(i)*xiray(i); 
            z = z + sstep*cray(i)*zetaray(i); 
            
            % add the new rray, zray values
            rray(i+1)=r; % the values you just calculated
            zray(i+1)=z;
            
            % Get new soundspeed for new location
            [cray(i+1), cz(i+1), czz(i+1)]=getCVal(rray(i+1),zray(i+1),cVal);
            
            zetaray(i+1) = zetaray(i) - sstep/(cray(i)^2)*cz(i); 
            xiray(i+1) = xiray(i); 
                  
            thetavec(i+1) = atan(zetaray(i+1)/xiray(i+1)); 
            
            % Move Forward
            sray(i+1)=sray(i)+sstep;
            dtau=sstep/cray(i);
            tau(i+1)=tau(i)+dtau;
            
            i=i+1;
            
        end
        tt = dth + 3*(dth-1);
        for ss = 1:length(thetavec)
            calc_new(tt,ss) = zray(1,ss);
            calc_new(tt+1,ss) = thetavec(1,ss);
            calc_new(tt+2,ss) = tau(1,ss);
            calc_new(tt+3,ss) = rray(1,ss);
        end
        save('calc_new', 'calc_new')
       
    end
end
function raytrace(sstep,depmat,cmat,thetas,zsources,f0,debug,z_glider)
% 2.681 Env. Ocean Acoustics HW 3

% sstep     the step size along the ray
% depmat    the depth matrix indexed by range, with form [r1,dep1;r2,dep2...]
% cmat      the soundspeed profile, with form [z1 c1;z2 c2...]
% thetas    the set of beam angles you want to run the raytrace code at,
%           between -90 and 90 degrees
% zsources  the source depths you want for beams
% f0        frequency
% numsteps  number of steps to run
% debug     flag that, if true, pauses at each top/bottom interaction

%% Interpolate Data
% If you've read COA, you'll know that you will need spatial derivatives of
% the sound speed field. It behooves you to calculate those derivatives
% first, and access them as nececessary, instead of making an interpolation
% every time.

[cVal,dVal] = dataInterp(cmat,depmat);

% dc = [cmat(:,1) gradient(cmat(:,2),cmat(:,1))]; 
% [czVal,~] = dataInterp(dc,depmat);
% 
% d2c = [cmat(:,1) gradient(dc(:,2),dc(:,1))];
% [czzVal,~] = dataInterp(d2c,depmat);
% 
% dz = [depmat(:,1) gradient(depmat(:,2),depmat(:,1))];
% [dzVal,~] = dataInterp(dz,depmat);

%Code Below does not work - DO NOT USE
%czVal = gradient(cVal.Values(:));
%czzVal = gradient(czVal(:));
%dzVal = gradient(dVal.Values(:));

%% Initialize Figures

% Raytrace Figure
raytraceFig=figure;
minR = 0; maxR = 0; % initialize for plotting purposes
hold on
set(gca,'YDir','Reverse')
xlabel('range (km)')
ylabel('depth (m)')
cc = hsv(length(thetas)*length(zsources));

thetavec_out = zeros(length(thetas), z_glider);
sray_out = zeros(length(thetas), z_glider);
rray_out = zeros(length(thetas), z_glider);
zray_out = zeros(length(thetas), z_glider);
tau_out = zeros(length(thetas), z_glider);


%% Ray Tracing
A0 = 1;
% loop over all sources
for dzs=1:length(zsources)
    z0=zsources(dzs); % get source depth
    pp = 1;
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
        q=[0];
        p=[1/c0];
        A=[A0]; %ray amplitude
        cz=[cz0]; %first derivative of ssp
        czz=[czz0]; %second derivative of ssp
        J0=0;
        i=1; %step count
        rsign=1;
        
        
        %% Iterate over numsteps to calculate out the ray
        
        while zray<z_glider
            
            r=rray(i);%r from previous step
            z=zray(i);%z from previous step
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%% INSERT CODE HERE TO GET r, z OF RAY IN NEW STEP %%%%%%
            
            r = r + sstep*cray(i)*xiray(i); 
            z = z + sstep*cray(i)*zetaray(i); 
           
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %figure out if there is bottom or surface interaction,
            %calculate new values
            
            [interact,sray,zray,rray,zetaray,xiray,cray,tau,thetavec,c0,theta0,cz,czz]=getInteraction(r,z,...
                dVal,sray,zray,rray,zetaray,xiray,cray,theta0,cVal,tau,thetavec,c0,cz,czz);
            
            %If there is an interaction, add two to 'i' and get new p, q, A
            if interact
                
                %add 2 to i to account for the two points from interaction fctn
                [q,p,A]=getAmplitudes(cray,zray,rray,sray(i)-sray(i-1),xiray,theta0,A0,i,p,q,A,c0,cz,czz); %First addition
                [q,p,A]=getAmplitudes(cray,zray,rray,sray(i+1)-sray(i),xiray,theta0,A0,i+1,p,q,A,c0,cz,czz); %Second
                i=i+2; %get off of the surface
                if debug %Plot the interaction if debug
                    plot(rray,zray)
                    hold on
                    plot(depmat(1:1000,1),depmat(1:1000,2))
                    hold off
                    figure
                    plot(sray(:,1:100),A(:,1:100))
                    pause
                end
                continue
            end
            
            % add the new rray, zray values
            rray(i+1)=r; % the values you just calculated
            zray(i+1)=z;
            
            % Get new soundspeed for new location
            [cray(i+1), cz(i+1), czz(i+1)]=getCVal(rray(i+1),zray(i+1),cVal);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%% PUT YOUR CODE TO GET zetaray(i+1), xiray HERE %%%%%%%
            
            zetaray(i+1) = zetaray(i) - sstep/(cray(i)^2)*cz(i); 
            xiray(i+1) = xiray(i); 
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%% PUT YOUR CODE TO GET NEW thetavec(i+1) HERE  %%%%%%%%%%
                        
            thetavec(i+1) = atan(zetaray(i+1)/xiray(i+1)); 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Move Forward
            sray(i+1)=sray(i)+sstep;
            dtau=sstep/cray(i);
            tau(i+1)=tau(i)+dtau;
            
            % next, calculate the amplitudes,p and q:
            [q,p,A]=getAmplitudes(cray,zray,rray,sstep,xiray,theta0,A0,i,p,q,A,c0,cz,czz);
            
            i=i+1;
            
            if debug
                if rem(i/100,1)==0
                  disp('on step number ')
                 i;
                end
            end
            
        end
        
        %% Plot Results
        
        figure(raytraceFig)
        plot(20,500);
        hold on 
        plot(rray/1000,zray,'color',cc(length(thetas)*(dzs-1)+dth,:),'linewidth',2)
        
        %save('output', 'thetavec', 'sray', 'rray', 'zray')
        for qq = 1:length(thetavec)
            thetavec_out(pp,qq) = thetavec(1,qq);
            sray_out(pp,qq) = sray(1,qq);
            rray_out(pp,qq) = rray(1,qq);
            zray_out(pp,qq) = zray(1,qq);
            tau_out(pp,qq) = tau(1,qq);
        end
        pp = pp+1;
        % Store min and max R for plotting purposes
        if min(rray) < minR
            minR = floor(min(rray));
        end
        
        if max(rray) > maxR
            maxR = ceil(max(rray));
        end
    end
    save('output', 'thetavec_out', 'sray_out', 'rray_out', 'zray_out', 'tau_out')
end

% make input bathymetry file fit full number of steps for plotting purposes
rSpace = minR:100:maxR;
bSpace = interp1(depmat(:,1),depmat(:,2), rSpace, 'linear', 'extrap');

rSpace(end+1) = rSpace(end);
bSpace(end+1) = max(bSpace)+100;

rSpace(end+1) = min(rSpace);
bSpace(end+1) = bSpace(end);

rSpace(end+1) = rSpace(1);
bSpace(end+1) = bSpace(1);

figure(raytraceFig)
plot(depmat(:,1)/1000,depmat(:,2),'color',1/255*[194 178 128]);
fill(rSpace/1000, bSpace, 1/255*[194 178 128],'linestyle','none');
axis tight
ylim([0 max(bSpace)]);
xlim([0 10]);

end
%% Prepare workspace
clear variables;
clc

%% Load data files
start_time=1576176163.79;
directory='D:\Documents\Thesis_Research\Bigelow_12_12\processed_data\PLATYPUS\';
%variables=["ACOUSTICRB","SOURCEXY", "NAVDR"]; % MOOS data only
variables=["PF_RB","NAVDR","ACOUSTICDATA"]; % Particle filter data
for ii=1:length(variables)
    filename=directory + string(start_time) + '_' + variables(ii) + '.mat';
    load(filename);
end
directory_ADCP="D:\Documents\Thesis_Research\Bigelow_12_12\processed_data\";
load(directory_ADCP + 'AD2CPData.346.00002_1_Proc.mat');

PF_RB=PF_RB(sel,:);
NAVDR(1,:)=[];
NAVDR=NAVDR(sel,:);
NAVDR.NAV_SPEED=NAVDR.NAV_SPEED-.5;

for jj=1:length(NAVDR.NAV_HEADING)
    if NAVDR.NAV_HEADING(jj)>180
        NAVDR.NAV_HEADING(jj)=NAVDR.NAV_HEADING(jj)-360;
    elseif NAVDR.NAV_HEADING(jj)<-180
        NAVDR.NAV_HEADING(jj)=NAVDR.NAV_HEADING(jj)+360;
    end
end
%% Define origin of local coordinate system
lat=41.524590;
lon=-70.671840;

%% Define UTM structure
zone=utmzone(lat,lon);
utmstruct = defaultm('utm'); 
utmstruct.zone = zone;  
utmstruct.geoid = wgs84Ellipsoid;
utmstruct = defaultm(utmstruct);

%% Convert source GPS data to UTM
[e_gps,n_gps]=mfwdtran(utmstruct,Data.jetyak_Lat,Data.jetyak_Lon);

%% Convert origin to UTM
[x,y] = mfwdtran(utmstruct,lat,lon);

%% Convert GPS to local coordinates
e_gps=e_gps-x;
n_gps=n_gps-y;

%% Interpolate data to common time base
total_time=PF_RB.Time(end)-PF_RB.Time(1);
tenth_second=0:.1:total_time;
Interp.Time=tenth_second+PF_RB.Time(1);

Interp.e_gps=interp1(Data.jetyak_Time,e_gps,Interp.Time);
Interp.n_gps=interp1(Data.jetyak_Time,n_gps,Interp.Time);
Interp.NAV_HEADING=interp1(NAVDR.TIME,NAVDR.NAV_HEADING,Interp.Time);
Interp.NAV_SPEED=interp1(NAVDR.TIME,NAVDR.NAV_SPEED,Interp.Time);
Interp.Range=interp1(PF_RB.Time,PF_RB.Range,Interp.Time);
Interp.Bearing=interp1(PF_RB.Time,PF_RB.Bearing,Interp.Time);


%% Define constants

%freq=10; %Number of cycles per second

% data_full=zeros(27,steps,iterations); %matrix for post-processing data
% epsilon_full=zeros(iterations,steps);
% epsilon_bar=zeros(1,steps);
% epsilon_v_full=zeros(iterations,steps);
% epsilon_v_bar=zeros(1,steps);
% mu_full=zeros(6,steps,iterations);
% mu_bar=zeros(6,steps);

%[data,epsilon_v]=NCV_field(NAVDR, SOURCEXY, ACOUSTICRB, e_gps_interp, n_gps_interp); % MOOS data
[data,epsilon_v,P_minus_out,H_out,S_out,K_out,z_k_out,h_out,F_out]=NCV_C_field(Interp); %PF data
% data_full(:,:,ii)=data;
% epsilon_full(ii,:)=epsilon;
% epsilon_bar=epsilon_bar+epsilon;
% epsilon_v_full(ii,:)=epsilon_v;
% epsilon_v_bar=epsilon_v_bar+epsilon_v;
% mu_full(:,:,ii)=mu;
% mu_bar=mu_bar+mu;
% 
% 
% epsilon_bar=epsilon_bar/iterations;
% epsilon_v_bar=epsilon_v_bar/iterations;
% mu_bar=mu_bar/iterations;
% 
% 
% %% Plot results
% chi_sq_upper1=chi2inv(0.975,6*iterations)/iterations;
% chi_sq_upper2=chi2inv(0.975,4*iterations)/iterations;
% chi_sq_lower1=chi2inv(0.025,6*iterations)/iterations;
% chi_sq_lower2=chi2inv(0.025,4*iterations)/iterations;
% 
% figure(1)
% plot(epsilon_bar)
% yline(chi_sq_upper1,'--b');
% yline(chi_sq_lower1,'--b');
% 
% figure(2)
% plot(epsilon_v_bar)
% yline(chi_sq_upper2,'--b');
% yline(chi_sq_lower2,'--b');
% 
% figure
% plot(1:steps,data_full(17,:,1))
% title('Error Covariance Matrix Norm vs Filter Step')
% xlabel('Filter Step')
% ylabel('Error Covariance Matrix Norm')
% 
% figure
% plot(1:steps,data_full(18,:,1))
% %axis([0 steps -5 5])
% title('Filter Error vs Filter Step')
% xlabel('Filter Step')
% ylabel('Filter Error (m)')
% 
% figure
% plot(1:steps,data_full(19,:,1))
% %axis([0 steps -5 5])
% title('Filter Error vs Filter Step')
% xlabel('Filter Step')
% ylabel('Filter Error (degrees)')
% 
% figure
% plot(data_full(7,:,1),data_full(9,:,1),data_full(15,:,1),data_full(16,:,1))
% axis equal
% title('Acutal Position vs Filter Output Position')
% xlabel('x position (m)')
% ylabel('y position (m)')
% legend('Filter Output','Actual Position','location','southeast')
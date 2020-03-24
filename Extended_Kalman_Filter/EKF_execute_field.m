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
%NAVDR.NAV_SPEED=NAVDR.NAV_SPEED-.5;

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


%% Execute filter
%[data,epsilon_v]=NCV_field(NAVDR, SOURCEXY, ACOUSTICRB, e_gps_interp, n_gps_interp); % MOOS data
[data,epsilon_v,rho_bar,P_minus_out,S_out,K_out,z_k_out,F_out]=NCV_C_field(Interp); %PF data

%% Plot results
steps=length(Interp.Time);
epsilon_v_bar_upper=chi2inv(0.995,4*steps)/steps;
epsilon_v_bar_lower=chi2inv(0.005,4*steps)/steps;
epsilon_v(isnan(epsilon_v))=0;
epsilon_v_bar=(1/steps)*sum(epsilon_v);
mag=hypot(data(11,:),data(12,:));
dir=atan2d(data(12,:),data(11,:));

plot_dir=90-dir;
for kk=1:length(plot_dir)
    if plot_dir(kk) <0
        plot_dir(kk)=plot_dir(kk)+360;
    end
end

figure(1)
plot(1:steps,data(15,:))
title('Error Covariance Matrix Norm vs Filter Step')
xlabel('Filter Step')
ylabel('Error Covariance Matrix Norm')

figure(2)
plot(data(7,:),data(9,:),z_k_out(1,:),z_k_out(3,:),Interp.e_gps,Interp.n_gps)
axis equal
title('Filter Output Position')
xlabel('Easting position (m)')
ylabel('Northing position (m)')
legend('Filter Output','Input Measurement','Jetyak Position')

figure(3)
reset(gca)
yyaxis left
plot(1:steps,data(16,:),1:steps,data(18,:))
title('Innovation vs Filter Step')
xlabel('Filter Step')
ylabel('Range (m)')
yyaxis right
plot(1:steps,data(17,:),1:steps,data(19,:))
ylabel('Speed (m/s)')
legend('Easting Range','Northing Range','Easting STW','Northing STW')

figure(4)
plot(1:steps,mag)
title('Calculated Current Magnitude vs Filter Step')
xlabel('Filter Step')
ylabel('Current Magnitude')

figure(5)
plot(1:steps,plot_dir)
title('Calculated Current Direction vs Filter Step')
xlabel('Filter Step')
ylabel('Current Direction')

figure(6)
reset(gca);
quiver(data(7,:),data(9,:),data(11,:),data(12,:))
title('Calculated Current Field')
xlabel('Easting position (m)')
ylabel('Northing position (m)')
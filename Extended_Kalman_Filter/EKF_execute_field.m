%% Prepare workspace
clear variables;
close all
clc

%% Load data files
start_time=1576176163.79; % 12/12
% start_time=1582217328.16; % 2/20
directory='D:\Documents\Thesis_Research\Bigelow_12_12\processed_data\PLATYPUS\';
% directory='D:\Documents\Thesis_Research\Bigelow_2_20\processed_data\platypus\';
%variables=["ACOUSTICRB","SOURCEXY", "NAVDR"]; % MOOS data only
% variables=["PF_RB","NAVDR","ACOUSTICDATA"]; % Particle filter data
variables=["MLE_RB","NAVDR","PF_RB","ACOUSTICDATA"]; % MLE data
for ii=1:length(variables)
    filename=directory + string(start_time) + '_' + variables(ii) + '.mat';
    load(filename);
end
directory_ADCP="D:\Documents\Thesis_Research\Bigelow_12_12\processed_data\";
load(directory_ADCP + 'AD2CPData.346.00002_1_Proc.mat'); % 12/12
% directory_ADCP="D:\Documents\Thesis_Research\Bigelow_2_20\processed_data\"; % 2/20
% load(directory_ADCP + 'AD2CPData.051.00000_Proc.mat'); % 2/20

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
one_second=0:1:total_time;
Interp.Time=one_second+PF_RB.Time(1);

[~,ind]=unique(MLE_RB.time);
MLE_RB=MLE_RB(ind,:);

Interp.e_gps=interp1(Data.jetyak_Time,e_gps,Interp.Time);
Interp.n_gps=interp1(Data.jetyak_Time,n_gps,Interp.Time);
Interp.NAV_HEADING=interp1(NAVDR.TIME,NAVDR.NAV_HEADING,Interp.Time);
Interp.NAV_SPEED=interp1(NAVDR.TIME,NAVDR.NAV_SPEED,Interp.Time);
Interp.Range=interp1(MLE_RB.time,MLE_RB.range,Interp.Time);
Interp.Bearing=interp1(MLE_RB.time,MLE_RB.bearing,Interp.Time);
Interp.Range_PF=interp1(PF_RB.Time,PF_RB.Range,Interp.Time);
Interp.Bearing_PF=interp1(PF_RB.Time,PF_RB.Bearing,Interp.Time);


%% Execute filter
%[data,epsilon_v]=NCV_field(NAVDR, SOURCEXY, ACOUSTICRB, e_gps_interp, n_gps_interp); % MOOS data
% [data,epsilon_v,rho_bar,P_minus_out,S_out,K_out,z_k_out,F_out]=NCV_C_field(Interp); %PF data
[data,epsilon_v,rho_bar,P_minus_out,S_out,K_out,z_k_out,F_out]=NCV_C_field_bias_hybrid(Interp); %PF data

%% Plot results
steps=length(Interp.Time);
epsilon_v_bar_upper=chi2inv(0.995,4*steps)/steps;
epsilon_v_bar_lower=chi2inv(0.005,4*steps)/steps;
epsilon_v(isnan(epsilon_v))=0;
epsilon_v_bar=(1/steps)*sum(epsilon_v);
mag=hypot(data(12,:),data(13,:));
dir=atan2d(data(13,:),data(12,:));

e_pos_PF=zeros(length(Interp.Bearing_PF),1);
n_pos_PF=zeros(length(Interp.Bearing_PF),1);

for mm=1:length(Interp.Bearing_PF)
    acoustic_bearing_PF=90-(Interp.NAV_HEADING(mm)-Interp.Bearing_PF(mm));
    if acoustic_bearing_PF>180
        acoustic_bearing_PF=acoustic_bearing_PF-360;
    elseif acoustic_bearing_PF<-180
        acoustic_bearing_PF=acoustic_bearing_PF+360;
    end
    source_x=Interp.Range_PF(mm)*cosd(acoustic_bearing_PF);
    source_y=Interp.Range_PF(mm)*sind(acoustic_bearing_PF);

    e_pos_PF(mm)=Interp.e_gps(mm)-source_x;
    n_pos_PF(mm)=Interp.n_gps(mm)-source_y;
end

plot_dir=90-dir;
for kk=1:length(plot_dir)
    if plot_dir(kk) <0
        plot_dir(kk)=plot_dir(kk)+360;
    end
end

figure(1)
plot(1:steps,data(17,:))
title('Error Covariance Matrix Norm vs Filter Step')
xlabel('Filter Step')
ylabel('Error Covariance Matrix Norm')

figure(2)
plot(data(8,:),data(10,:),z_k_out(1,:),z_k_out(2,:),e_pos_PF,n_pos_PF,Interp.e_gps,Interp.n_gps)
axis equal
title('Filter Output Position')
xlabel('Easting position (m)')
ylabel('Northing position (m)')
legend('Filter Output','Input Measurement','PF Position','Jetyak Position')

figure(3)
reset(gca)
yyaxis left
plot(1:steps,data(18,:),1:steps,data(19,:))
title('Innovation vs Filter Step')
xlabel('Filter Step')
ylabel('Range (m)')
yyaxis right
plot(1:steps,data(20,:),1:steps,data(21,:))
ylabel('Speed (m/s)')
legend('Easting Range','Northing Range','STW','Heading')

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
quiver(data(8,:),data(10,:),data(12,:),data(13,:))
title('Calculated Current Field')
xlabel('Easting position (m)')
ylabel('Northing position (m)')
%% Plot ADCP current
if strcmpi(Config.Burst_BottomTrack, 'True')
    dataModeWord = 'Comp';
elseif strcmpi( Config.Plan_BurstEnabled, 'True' )
    dataModeWord = 'Burst';
elseif strcmpi( Config.Plan_AverageEnabled, 'True' )
    dataModeWord = 'Average';
end

if strcmpi(Config.Burst_BottomTrack, 'True')
    planModeWord='BottomTrack';
else
    planModeWord=dataModeWord;
end

avg_VelEast=zeros(length(Data.( [ dataModeWord '_VelEast' ] )),1);
avg_VelNorth=zeros(length(Data.( [ dataModeWord '_VelEast' ] )),1);
for ii=1:length(Data.( [ dataModeWord '_VelEast' ] ))
    avg_VelEast(ii)=mean(Data.( [ dataModeWord '_VelEast' ] )(ii,:));
    avg_VelNorth(ii)=mean(Data.( [ dataModeWord '_VelNorth' ] )(ii,:));
end



Interp.avg_VelEast=interp1(Data.( [ planModeWord '_Time_UNIX' ] ),avg_VelEast,Interp.Time);
Interp.avg_VelNorth=interp1(Data.( [ planModeWord '_Time_UNIX' ] ),avg_VelNorth,Interp.Time);

mag_avg=hypot(Interp.avg_VelEast,Interp.avg_VelNorth);
dir_avg=atan2d(Interp.avg_VelNorth,Interp.avg_VelEast);

% figure(7)
% yyaxis left
% plot(Interp.Time,mag_avg)
% yyaxis right
% plot(Interp.Time,Interp.NAV_HEADING)
% 
% figure(8)
% yyaxis left
% plot(Interp.Time,dir_avg)
% yyaxis right
% plot(Interp.Time,Interp.NAV_HEADING)

figure(9)
reset(gca);
quiver(Interp.e_gps,Interp.n_gps,Interp.avg_VelEast,Interp.avg_VelNorth,2)

scale_factor=10;
figure(10)
reset(gca)

quiver(data(8,:),data(10,:),scale_factor*data(12,:),scale_factor*data(13,:),'AutoScale','off')
hold on
quiver(Interp.e_gps,Interp.n_gps,scale_factor*Interp.avg_VelEast,scale_factor*Interp.avg_VelNorth,'AutoScale','off')
hold on
quiver(-100,-100,scale_factor*0,scale_factor*1,'AutoScale','off')
axis equal
legend('Calculated current','Measured current','1 m/s')
xlabel('Easting position (m)')
ylabel('Northing position (m)')
%% Prepare workspace
clear variables
clc

%% Load data files
% % 10/30/19
% directory=["D:\Documents\Thesis_Research\Bigelow_10_30\raw_data\"];
% directory1='D:\Documents\Thesis_Research\Bigelow_10_30\processed_data\';
% files_ADCP=["AD2CPData.303.00000"];
% files_jetyak_GPS=["00000006.log-7209129"];
% % 12/12/19
% directory=["D:\Documents\Thesis_Research\Bigelow_12_12\raw_data\"];
% directory1='D:\Documents\Thesis_Research\Bigelow_12_12\processed_data\';
% files_ADCP=["AD2CPData.346.00002_1"];
% files_jetyak_GPS=["00000003.BIN-11044046"];
% 2/20/20
directory=["D:\Documents\Thesis_Research\Bigelow_2_20\raw_data\"];
directory1='D:\Documents\Thesis_Research\Bigelow_2_20\processed_data\';
files_ADCP=["AD2CPData.051.00000"];
files_jetyak_GPS=["00000013.log-8869768"];
%% Parse GPS data
filename=directory(1) + files_jetyak_GPS(1)+'.mat';
load(filename);
ind =1;
while GPS(ind,3)<3
    ind=ind+1;
end
GPS_trimmed=GPS(ind:end,:);
jetyak_time=zeros(length(GPS_trimmed),1);
for ii=1:length(GPS_trimmed)
    jetyak_time(ii)=315964782+((GPS_trimmed(ii,5)*7*86400)+(GPS_trimmed(ii,4)*.001));
end

%% Initialize constants
twoZs=0; % 0 if false, 1 if true
% Correct for magnetic declination
declination_degrees=14;
declination_minutes=23;
declination_direction='west'; %'west' or 'east'
declination_correction=declination_degrees+(declination_minutes/60);
if strcmpi( declination_direction, 'east' )
	declination_correction=declination_correction*-1;
end

%% Process files

filename=directory(1) + files_ADCP(1)+'.mat';
load(filename);
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

% Define beam angles
Config.BeamCfg1_theta=25;
Config.BeamCfg2_theta=25;
Config.BeamCfg3_theta=25;
Config.BeamCfg4_theta=25;


Config.BeamCfg1_phi=0;
Config.BeamCfg2_phi=-90;
Config.BeamCfg3_phi=180;
Config.BeamCfg4_phi=90;

Config.BeamCfg5_theta=0;
Config.BeamCfg5_phi=0;

cfl= 13; %number of points in median filter for cleaning Bottom Track data

%% Interpolate Burst/Average data to Bottom Track time base
if strcmpi(Config.Burst_BottomTrack, 'True')
    for ii=1:Data.Burst_NBeams(1,1)
        Data.( [ 'BottomTrack_VelBeam' num2str( ii ) '_Clean'] )= clean0(Data.( [ 'BottomTrack_VelBeam' num2str( ii )] ), cfl, 1); 
        Data.( [ 'Burst_VelBeam' num2str( ii ) '_Interp' ] )=interp1(Data.Burst_Time, Data.( [ 'Burst_VelBeam' num2str( ii )] ), Data.BottomTrack_Time);
        Data.( [ 'Comp_VelBeam' num2str( ii )] )=Data.( [ 'Burst_VelBeam' num2str( ii ) '_Interp' ] )-Data.( [ 'BottomTrack_VelBeam' num2str( ii ) '_Clean'] );
    end
end

%% Perform compass calibration
% calibrate compass to adjust for hard iron of motor. can adjust caibration
Par         = CircleFitByPratt([Data.( [ planModeWord '_Magnetometer' ] )(:,1),...
              Data.( [ planModeWord '_Magnetometer' ] )(:,2)]); % fit magnetometer data
ang         = 0:0.01:2*pi; 
xp          = Par(3)*cos(ang); yp=Par(3)*sin(ang);

% calibrate compass based on magentometer offset in Par 
Data      = CalibrateCompass8_2017(Data, Par, planModeWord); % calibrate compass for hard iron correction from Jetyak motor

% % plot calibrated magnetometer data
% figure(1);clf
% scatter(Data.( [ planModeWord '_Magnetometer' ] )(:,1), Data.( [ planModeWord '_Magnetometer' ] )(:,2), 'k'); hold on
% plot(Par(1)+xp,Par(2)+yp, 'r','LineWidth',2);
% scatter(Data.( [ planModeWord '_Magnetometer_Cal' ] )(:,1), Data.( [ planModeWord '_Magnetometer_Cal' ] )(:,2),'b.')
% legend('Original', 'Fit', 'New','Location','NorthWest');
% title(['HxHy = [',num2str(-Par(1),'%5.2f'),' ',num2str(-Par(2),'%5.2f'),']']);
% axis equal; grid; xlabel('Magnetometer X'); ylabel('Magnetometer Y');
% %scatter(Data.Burst_MagnetometerX(ind), Data.Burst_MagnetometerY(ind), 'r'); hold on
%% Convert from beam to XYZ and ENU coordinate systems
[ Data, Config, T_beam2xyz ] = signatureAD2CP_beam2xyz_enu( Data, Config, dataModeWord, planModeWord, twoZs, declination_correction);

%% Convert to UNIX time
Data.( [ planModeWord '_Time_UNIX' ] ) = (Data.( [ planModeWord '_Time' ] )-719529)*(24*60*60);

%% Calculate magnitude and direction and plot vs time
mag=sqrt(Data.( [ dataModeWord '_VelEast' ] ).^2+Data.( [ dataModeWord '_VelNorth' ] ).^2);
dir=atan2d(Data.( [ dataModeWord '_VelNorth' ] ),Data.( [ dataModeWord '_VelEast' ] ));

avg_mag=zeros(length(Data.( [ dataModeWord '_VelEast' ] )),1);
avg_dir=zeros(length(Data.( [ dataModeWord '_VelEast' ] )),1);
for ii=1:length(Data.( [ dataModeWord '_VelEast' ] ))
    avg_mag(ii)=mean(mag(ii,:));
    avg_dir(ii)=mean(dir(ii,:));
end

avg_VelEast=zeros(length(Data.( [ dataModeWord '_VelEast' ] )),1);
avg_VelNorth=zeros(length(Data.( [ dataModeWord '_VelEast' ] )),1);
for ii=1:length(Data.( [ dataModeWord '_VelEast' ] ))
    avg_VelEast(ii)=mean(Data.( [ dataModeWord '_VelEast' ] )(ii,:));
    avg_VelNorth(ii)=mean(Data.( [ dataModeWord '_VelNorth' ] )(ii,:));
end

mag_avg=hypot(avg_VelEast,avg_VelNorth);
dir_avg=atan2d(avg_VelNorth,avg_VelEast);

% figure
% yyaxis left
% plot(Data.( [ planModeWord '_Time_UNIX' ] ),mag(:,3))
% yyaxis right
% plot(Data.( [ planModeWord '_Time_UNIX' ] ),Data.( [ planModeWord '_Heading_Cal' ] ))
% 
% figure
% yyaxis left
% plot(Data.( [ planModeWord '_Time_UNIX' ] ),dir(:,3))
% yyaxis right
% plot(Data.( [ planModeWord '_Time_UNIX' ] ),Data.( [ planModeWord '_Heading_Cal' ] ))
% 
% figure
% yyaxis left
% plot(Data.( [ planModeWord '_Time_UNIX' ] ),mag_avg)
% yyaxis right
% plot(Data.( [ planModeWord '_Time_UNIX' ] ),Data.( [ planModeWord '_Heading' ] ))
% 
% figure
% yyaxis left
% plot(Data.( [ planModeWord '_Time_UNIX' ] ),dir_avg)
% yyaxis right
% plot(Data.( [ planModeWord '_Time_UNIX' ] ),Data.( [ planModeWord '_Heading' ] ))

figure
yyaxis left
plot(Data.( [ planModeWord '_Time_UNIX' ] ),avg_mag)
yyaxis right
plot(Data.( [ planModeWord '_Time_UNIX' ] ),Data.( [ planModeWord '_Heading_Cal' ] ))

figure
yyaxis left
plot(Data.( [ planModeWord '_Time_UNIX' ] ),avg_dir)
yyaxis right
plot(Data.( [ planModeWord '_Time_UNIX' ] ),Data.( [ planModeWord '_Heading_Cal' ] ))

%% Interpolate GPS data to ADCP time base
Data.jetyak_Time=jetyak_time;
Data.jetyak_Lat=GPS_trimmed(:,8);
Data.jetyak_Lon=GPS_trimmed(:,9);
Data.jetyak_Lat_Interp=interp1(jetyak_time, GPS_trimmed(:,8),Data.( [ planModeWord '_Time_UNIX' ] ));
Data.jetyak_Lon_Interp=interp1(jetyak_time, GPS_trimmed(:,9),Data.( [ planModeWord '_Time_UNIX' ] ));
 %% Plot quiver
reset(gca);
quiver(Data.jetyak_Lon_Interp(2300:21120),Data.jetyak_Lat_Interp(2300:21120),avg_VelEast(2300:21120),avg_VelNorth(2300:21120),2)
%% Save updated Data struct
save(directory1 + files_ADCP(1)+ '_Proc.mat','Config','Data','Descriptions','Units','-v7.3')
% This file processes S1000 AD2CP data from the Jetyak and combines it with
% GPS data. 
% Katie Samuelson
% Code written 8/15/2014 (former name: ProcessADCPwithGPS(_WK))
% Wouter Kranenburg; may/june 2017;
% new name: adcp_2ENU


clear 
clc
close all
addpath('Functions')

if 0
% load ADCP data 
daynr       = 212;
files2load  = dir(['Data.mat\*',num2str(daynr),'*.mat']);                            % .mat files to load from yearday 236
DataAllSep  = [];                                                           % concatenate the structures from each data file
fnms        = [];
imax        = length(files2load);
for i = 1:imax
    loadfile    = files2load(i);
    load([loadfile.folder,'\',loadfile.name]);
    DataAllSep  = [DataAllSep, Data];
    f           = fieldnames(Data);                                         % field names
    fnms        = [fnms, f];
end
if imax>1                                                                   %action needed to concatenate 
    C       = reshape(fieldnames(Data), 1, []); 
%     DataAll = struct(C{:});                                                 % new data structure that will concatenate each field 
    for j   = 1:length(f)                                                   % concatenate each field within the structures
        fieldtouse = char(f(j)); 
        if daynr==206 %now rows for 206
            DataAll.(fieldtouse) = cat(1,DataAllSep.(fieldtouse));
        else
%        if j==1||j==2||j==40||j==43||j==44||j==45||j==74||j==77||j==78||j==79||j==118 % original line KS
%        if j==1||j==2||j==44||j==48||j==49||j==50||j==88||j==89
         if j==1||j==2||j==48||j==49||j==86||j==87
             disp([num2str(j),'-->',char(f(j))]); %write on screen
             DataAll.(fieldtouse) = cat(2,DataAllSep.(fieldtouse)); %these fields are rows iso columns; Mostly Time Related
         else
             DataAll.(fieldtouse) = cat(1,DataAllSep.(fieldtouse));
         end 
        end
        f2 = fieldnames(DataAll);
    end
    [fnms, f2]                                                             %check fieldnames; apparently order of fields in new struct different
    Data = DataAll;
end

tekst1 = datestr(Data.Burst_MatlabTimeStamp(1),'yyyymmddHHMM');
tekst2 = datestr(Data.Burst_MatlabTimeStamp(end),'yyyymmddHHMM');
tekst = ['Data Time Range = ',char(tekst1),' to ',char(tekst2)];
disp(tekst)
end


% load GPS data
% load('Data\20140824a_try_01_varbs'); 
load Data.137.00002.mat
%% Compass calibration
% calibrate compass to adjust for hard iron of motor. can adjust caibration
Par         = CircleFitByPratt([Data.Burst_MagnetometerX,...
              Data.Burst_MagnetometerY]); % fit magnetometer data
ang         = 0:0.01:2*pi; 
xp          = Par(3)*cos(ang); yp=Par(3)*sin(ang);

% calibrate compass based on magentometer offset in Par 
DataNew      = CalibrateCompass8_2017(Data, Par); % calibrate compass for hard iron correction from Jetyak motor
Heading_ADCP = DataNew.Burst_Heading+16.1-90; % 16.1 accounts for magnetic declination

% plot calibrated magnetometer data
figure(1);clf
scatter(Data.Burst_MagnetometerX, Data.Burst_MagnetometerY, 'k'); hold on
plot(Par(1)+xp,Par(2)+yp, 'r');
scatter(DataNew.Burst_MagnetometerX, DataNew.Burst_MagnetometerY,'b')
legend('Original', 'Fit', 'New');
title(['HxHy = [',num2str(-Par(1),'%5.2f'),' ',num2str(-Par(2),'%5.2f'),']']);
axis equal; grid; xlabel('Magnetometer X'); ylabel('Magnetometer Y');
%saveas(gcf,['Data.calibrated/CompassCal',num2str(daynr),'.png'],'png');

%% Time Check (and assignment for day 206) & Calibration with ADCP-up excluded
%the above doesn't seem to work? Yes. Question is: is it good enough? (WK)
% FIND WHERE in timeseries IS THE BIG MAGNETOMETER DEVIATION
index   = 1:length(Data.Burst_TimeStamp);
ind     = index(abs(Data.Burst_MagnetometerX)>300);

if daynr==206
    Data.Burst_HostTimeMatlab  = Data.Burst_MatlabTimeStamp';
    Data.Burst_HostTime        = Data.Burst_TimeStamp' ;
    Data.BurstBT_HostTimeMatlab= Data.BurstBT_MatlabTimeStamp';
    Data.BurstBT_HostTime      = Data.BurstBT_TimeStamp'  ; 
    Data.IBurst_HostTimeMatlab = Data.IBurst_MatlabTimeStamp';
    Data.IBurst_HostTime       = Data.IBurst_TimeStamp' ;     
end

clear B BT IB
B(1,:)=   (Data.Burst_HostTimeMatlab   - Data.Burst_MatlabTimeStamp')*24*3600;
B(2,:)=    Data.Burst_HostTime         - Data.Burst_TimeStamp' ;
BT(1,:)=  (Data.BurstBT_HostTimeMatlab - Data.BurstBT_MatlabTimeStamp')*24*3600;
BT(2,:)=   Data.BurstBT_HostTime       - Data.BurstBT_TimeStamp'  ; 
IB(1,:)=  (Data.IBurst_HostTimeMatlab  - Data.IBurst_MatlabTimeStamp')*24*3600;
IB(2,:)=   Data.IBurst_HostTime        - Data.IBurst_TimeStamp' ;     
figure
subplot(2,2,1); hh = plot(1:size(B,2),B); title('diff B'); set(hh(1),'LineWidth',1.8); set(hh(2),'LineStyle','--'); grid on; ylabel('[s]');
subplot(2,2,2); hh = plot(1:size(BT,2),BT); title('diff BT');set(hh(1),'LineWidth',1.8); set(hh(2),'LineStyle','--'); grid on; ylabel('[s]'); 
subplot(2,2,3); hh =plot(1:size(IB,2),IB); title('diff IB');set(hh(1),'LineWidth',1.8); set(hh(2),'LineStyle','--'); grid on; ylabel('[s]');
annotation('textbox',[0.6 0.2 0.35 0.1],'LineStyle','none','String',['HostTime minus TimeStamp for Burst,BurstBT,IBurst'],'FontSize',10);
saveas(gcf,['Data.calibrated/TimeStamps',num2str(daynr),'.png'],'png');

xdata = [daynr-1+datenum('20170101','yyyymmdd'):1/24:daynr+datenum('20170101','yyyymmdd')];
t1    = Data.Burst_HostTimeMatlab(1); 

figure
subplot(2,1,1)
plot(Data.Burst_HostTimeMatlab -t1    ,Data.Burst_MagnetometerX,'k.'); hold on
plot(Data.Burst_HostTimeMatlab(ind)-t1,Data.Burst_MagnetometerX(ind),'r.'); hold on
ylabel('MagnetometerX');

subplot(2,1,2)
plot(Data.Burst_HostTimeMatlab-t1     ,Data.Burst_MagnetometerY,'k.'); hold on
plot(Data.Burst_HostTimeMatlab(ind)-t1,Data.Burst_MagnetometerY(ind),'r.'); hold on
ylabel('MagnetometerY');
saveas(gcf,['Data.calibrated/CompassCal',num2str(daynr),'_errors.png'],'png');

figure(1)
scatter(Data.Burst_MagnetometerX(ind), Data.Burst_MagnetometerY(ind), 'r'); hold on
saveas(gcf,['Data.calibrated/CompassCal',num2str(daynr),'.png'],'png');

%% NEW Compass calibration
index   = 1:length(Data.Burst_TimeStamp);
tbu     = index(abs(Data.Burst_MagnetometerX)<=300); %so: which to INCLUDE (to be used)

% calibrate compass to adjust for hard iron of motor. can adjust caibration
Par         = CircleFitByPratt([Data.Burst_MagnetometerX(tbu),...
              Data.Burst_MagnetometerY(tbu)]); % fit magnetometer data
ang         = 0:0.01:2*pi; 
xp          = Par(3)*cos(ang); yp=Par(3)*sin(ang);

% calibrate compass based on magentometer offset in Par 
DataNew      = CalibrateCompass8_2017(Data, Par); % calibrate compass for hard iron correction from Jetyak motor
Heading_ADCP = DataNew.Burst_Heading+16.1-90; % 16.1 accounts for magnetic declination

% plot calibrated magnetometer data
figure(1);clf
scatter(Data.Burst_MagnetometerX, Data.Burst_MagnetometerY, 'k'); hold on
plot(Par(1)+xp,Par(2)+yp, 'r','LineWidth',2);
scatter(DataNew.Burst_MagnetometerX, DataNew.Burst_MagnetometerY,'b.')
legend('Original', 'Fit', 'New','Location','NorthWest');
title(['HxHy = [',num2str(-Par(1),'%5.2f'),' ',num2str(-Par(2),'%5.2f'),']']);
axis equal; grid; xlabel('Magnetometer X'); ylabel('Magnetometer Y');
scatter(Data.Burst_MagnetometerX(ind), Data.Burst_MagnetometerY(ind), 'r'); hold on
saveas(gcf,['Data.calibrated/CompassCal',num2str(daynr),'_new.png'],'png');



%% Rotate from BEAM to XYZ and ENU
%Dataxyz, Configxyz, T_beam2xyz] = signatureAD2CPr2226_beam2xyz_enu8_2014(DataNew, Config, 'bt');
[Dataxyz, Configxyz, T_beam2xyz] = signatureAD2CPr2226_beam2xyz_enu8_2017(DataNew, Config, 'bt');
[Dataxyz, Configxyz, T_beam2xyz] = signatureAD2CPr2226_beam2xyz_enu8_2017(Dataxyz, Configxyz, 'burst');

% clean bottom track data
clear B BT IB
cfl         = 13;
BT.Vx       = clean0(Dataxyz.BurstBT_VelX    , cfl, 1); 
BT.Vy       = clean0(Dataxyz.BurstBT_VelY    , cfl, 1);
BT.Vz       = clean0(Dataxyz.BurstBT_VelZ    , cfl, 1);
BT.VEast    = clean0(Dataxyz.BurstBT_VelEast , cfl, 1);
BT.VNorth   = clean0(Dataxyz.BurstBT_VelNorth, cfl, 1);
BT.VUp      = clean0(Dataxyz.BurstBT_VelUp   , cfl, 1);

%% Interpolation time axis

z = [Config.burst_blankingDistance:Config.burst_cellSize:...
    Config.burst_nCells*Config.burst_cellSize-Config.burst_blankingDistance]; % AD2CP cell locations

% velocity profile excluding Jetyak motion % need to interpolate between the bottom track to be at the burst data times. 
% if daynr~=206% original lines: if ~206
    BurstTimeStamps = Dataxyz.Burst_HostTime; 
    BTTimeStamps    = Dataxyz.BurstBT_HostTime;
% else
%     BurstTimeStamps = Dataxyz.Burst_TimeStamp; 
%     BTTimeStamps    = Dataxyz.BurstBT_TimeStamp; %This one seems not OK
% end    

for i= 1:length(z)
VelE(i,1:length(BTTimeStamps)) = interp1(BurstTimeStamps, Dataxyz.Burst_VelEast(:,i) , BTTimeStamps)';
VelN(i,1:length(BTTimeStamps)) = interp1(BurstTimeStamps, Dataxyz.Burst_VelNorth(:,i), BTTimeStamps)';
VelU(i,1:length(BTTimeStamps)) = interp1(BurstTimeStamps, Dataxyz.Burst_VelUp(:,i)   , BTTimeStamps)';
end

VEastDiff   = VelE - repmat(BT.VEast,  [1 20])';
VNorthDiff  = VelN - repmat(BT.VNorth, [1 20])';
VUpDiff     = VelU - repmat(BT.VUp,    [1 20])';

% average upper water column velocity
VEastDiffMeanUpper  = mean(VEastDiff (2:5, :));
VNorthDiffMeanUpper = mean(VNorthDiff(2:5, :));


%% Velocity profiles; NB: filtering/smoothing?
MatlabTimeStampCorrected    = Data.Burst_MatlabTimeStamp;
day1                        = MatlabTimeStampCorrected(1);
filt                        = ones(1,4)./4;
fVEastDiff                  = imfilter(VEastDiff,filt);
fVNorthDiff                 = imfilter(VNorthDiff,filt);
fVUpDiff                    = imfilter(VUpDiff,filt);

%figure Vel EastNorthUp 
figure(5);clf
subplot(311)
imagesc((MatlabTimeStampCorrected-day1)*24*60, z, fVEastDiff); 
c = colorbar('peer', gca); set(get(c, 'ylabel'), 'String', 'Velocity (m/s)', 'Rotation', -90);
caxis([-1 1]*.5);
title('Velocity East'); xlabel('Time (min)'); ylabel('Depth from Surface (m)')

subplot(312)
imagesc((MatlabTimeStampCorrected-day1)*24*60, z, fVNorthDiff)
c = colorbar('peer', gca); set(get(c, 'ylabel'), 'String', 'Velocity (m/s)', 'Rotation', -90);
caxis([-.75 .5]);
title('Velocity North'); xlabel('Time (min)'); ylabel('Depth from Surface (m)')

subplot(313)
imagesc((MatlabTimeStampCorrected-day1)*24*60, z, fVUpDiff)
c = colorbar('peer', gca); set(get(c, 'ylabel'), 'String', 'Velocity (m/s)', 'Rotation', -90);
caxis([-0.25 0.25]);
title('Velocity Up'); xlabel('Time (min)'); ylabel('Depth from Surface (m)')

% yaxis([0 8]);yzoom
% xaxis([14.4 15.5]);xzoom
for i = 1:3
    subplot(3,1,i)
    ylim([0 8]);
%     xlim([100 150])
end

% save (turned off, takes time)
save(['Data.calibrated/Data.',num2str(daynr),'.cal.mat'])


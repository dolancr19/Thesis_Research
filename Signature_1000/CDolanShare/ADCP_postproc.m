% name:   checkBT_GPS
% author: W.M.Kranenburg
% date:   April/May/June 2017
% function: 
% - check BT with GPS

% pre-requisite
% - ADCP: data2ENU
% - GPS: nmea2mat
% ------------------------------

for yrday = [108, 109,117,137,138,205,206,209,212];

clearvars -except yrday
close all
clc
addpath('Functions');
drpbx     = 'c:\users\wmkra\Dropbox\';
addpath([drpbx,'mrocky']); 
addpath([drpbx,'mwouter']);
fldrgps   = [drpbx,'north_river\north_river_2017\Sig1000data\GPS\'];


%% Initialization: load GPS and ADCP
% yrday    = 137;                                                             %nb: this is the 'ADCP-yearday', actually 1 too high

date     = yrday+datenum('20170101','yyyymmdd')-1;
% name GPS-file tbu
if date<datenum('20170601','yyyymmdd');
    GPSfiles = dir([fldrgps,'GPS-',datestr(date,'yyyymmdd'),'-*.mat']);     %JetyakGPS
elseif date==datenum('20170724','yyyymmdd') || date==datenum('20170728','yyyymmdd')
    GPSfiles = dir([fldrgps,'NR',datestr(date,'mmddyy'),'*.mat']);         %GPSzodiac
elseif date==datenum('20170725','yyyymmdd') || date==datenum('20170731','yyyymmdd')
    GPSfiles = dir([fldrgps,datestr(date,'yyyymmdd'),'*.mat']);            %GPSbase
end
    
load(['Data.calibrated/Data.',num2str(yrday),'.cal.mat']);      %ADCP
GPS      = load([GPSfiles(1).folder,'\',GPSfiles(1).name]);     %GPS

% 'corrections' wrt GPS (or to get it in same format/pars as used earlier)
if date<datenum('20170601','yyyymmdd') %all Jetyak survey
    disp('using JetYak GPS info');
    % time axes for GPS: 1 correction and 1 extra axis (timeshort, speed and truecourse on shorter axis)
    ind   = find(diff(GPS.time)==0);
    if ~isempty(ind); GPS.time(ind+1) = nan; end

    if GPS.Ldif == 2
        GPS.timeshort = GPS.time(2:end-1);
    elseif GPS.Ldif == 4
        GPS.timeshort = GPS.time(3:end-2);
    elseif GPS.Ldif == 5
        GPS.timeshort = GPS.time(4:end-2);    
    end
else 
    if GPSfiles.name(10:12)=='zod' %all zodiac-gps-files
        disp('using zodiac GPS info');
        disp('switching GPS back to local time');
        GPS.dn  = GPS.dn-4/24;
%         depth   = interp(GPS.dn(2:end),GPS.depth,dn); %(dn could be dn-adcp
%         depGPS  = depth;
    elseif GPSfiles.name(10:13)=='base' %all used base-gps-files
        disp('using non-zodiac GPS info');
        disp('time already set to local in nmea2mat');
        GPS.dn = GPS.time;
        indnonnanlat = find(~isnan(GPS.lat)); %remove nans
        GPS.dn  = GPS.dn (indnonnanlat);
        GPS.lat = GPS.lat(indnonnanlat);
        GPS.lon = GPS.lon(indnonnanlat);
    end
    GPS.time = GPS.dn;
    GPS.tstart = GPS.time(1);
    GPS.tend = GPS.time(end);
    plot(GPS.dn,GPS.lon,'r.')
end


% time axes from ADCP
diffUTCEST  = 4; %in hrs
timeBT      = Data.BurstBT_HostTimeMatlab - diffUTCEST/24;
% timemin     = (MatlabTimeStampCorrected(1:end)-day1)*24*60; %in min from start ADCP  %sometimes (1:end-1)! (length DistBed)
timemin     = (timeBT-timeBT(1))*24*60;



%% Procedures
PROC01_CheckBT                                                             % Check BT with GPS
PROC01_CheckBT_figs                                                        % Further improvement BT possible/needed --> change e.g. heading-cor in calibration
PROC02_BLandmask                                                           % STEP 2: Make mask using BL from BT

dep = interp1(timeBT,DistBed2,GPS.time);                                    % dep used in Indexing for picture


Result_Indexing = ['Figures\',num2str(yrday),'\04_Indexing\TABELlong.mat'];
if exist(Result_Indexing)~=2;
    PROC03_Indexing                                                             % split up in lateral transect/diag/etc; small manual_correction
    PROC03_Indexing_FigsTab
    PROC04_IndexingLong                                                         % do the same on long time axis (improvement possible)
    PROC04_IndexingLong_FigsTab
else
    load(Result_Indexing); %results from earlier run on Indexing
    message = 'TABELlong does exist'
end


% make sure that below bottom is excluded in determination max direction etc. (should actually be ok, but to be sure)
uE2            = zeros(size(uE));
uE2(mask==1)   = uE2(mask==1)+uE(mask==1);
uN2            = zeros(size(uE));
uN2(mask==1)   = uN2(mask==1)+uN(mask==1);


% PROC05_Arrows                     %figs only
PROC06_Rotation                     %proc (theta&rot)
PROC06_Rotation_Figs              %figs only
% PROC07_VelocityFigs               %figs only
% PROC07b_VelocityFigs_rotfixed     %figs only
% PROC08_FigsGathered               %figs only

% PROC09_ImproveIndex               %proc to better split Up&Down, and make yNR(ind) monotone --> results added to L(.tbu)
% PROC10_FigsGathered               %figs only, just new ones for the improved cases

PROC11_SaveStuff

end









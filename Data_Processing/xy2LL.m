%% Prepare workspace
clear variables
clc

%% Define start time
%start_time=1574798368.75; %11-26
start_time=1576176163.79; %12-12 Platypus 1
%start_time=1576179399.92; %12-12 Platypus 2
%start_time=1576180456.34; %12-12 Platypus 3
%start_time=1576176663.37; %12-12 Quokka 1, no submerged time
%start_time=1576177784.93; %12-12 Quokka 2
%start_time=1576178198.67; %12-12 Quokka 3

%% Load data files
directory="D:\Documents\Thesis_Research\Bigelow_12_12\processed_data\PLATYPUS\";
%directory="D:\Documents\Thesis_Research\Bigelow_12_12\processed_data\QUOKKA\";
directory_ADCP="D:\Documents\Thesis_Research\Bigelow_12_12\processed_data\";
load(directory + string(start_time) + '_SOURCEXY.mat');
load(directory + string(start_time) + '_NAVDEPTH.mat');
load(directory + string(start_time) + '_NAVLL.mat');
load(directory_ADCP + 'AD2CPData.346.00002_1_Proc.mat');

%% Define origin of local coordinate system
lat=41.524590;
lon=-70.671840;

%% Define UTM structure
zone=utmzone(lat,lon);
utmstruct = defaultm('utm'); 
utmstruct.zone = zone;  
utmstruct.geoid = wgs84Ellipsoid;
utmstruct = defaultm(utmstruct);

%% Convert origin to UTM
[x,y] = mfwdtran(utmstruct,lat,lon);

%% Convert source GPS data to UTM
[x_gps,y_gps]=mfwdtran(utmstruct,Data.jetyak_Lat,Data.jetyak_Lon);
x_gps_interp=interp1(Data.jetyak_Time,x_gps,SOURCEXY.TIME);
y_gps_interp=interp1(Data.jetyak_Time,y_gps,SOURCEXY.TIME);
depth_interp=interp1(NAVDEPTH.TIME,NAVDEPTH.NAV_DEPTH,SOURCEXY.TIME);
depth_interp(isnan(depth_interp))=0;
navLat_interp=interp1(NAVLL.TIME,NAVLL.NAV_LAT,SOURCEXY.TIME);
navLon_interp=interp1(NAVLL.TIME,NAVLL.NAV_LONG,SOURCEXY.TIME);

%% Convert XY data to UTM

data_utm=[x_gps_interp-SOURCEXY.ACOUSTIC_SOURCE_X y_gps_interp-SOURCEXY.ACOUSTIC_SOURCE_Y];
data_jetyak=[x_gps_interp y_gps_interp];
%% Convert UTM to LL
[lat1,lon1] = minvtran(utmstruct,data_utm(:,1),data_utm(:,2));
[jetyakLat_interp,jetyakLon_interp]=minvtran(utmstruct,data_jetyak(:,1),data_jetyak(:,2));
%data_LL=[lat1 lon1];

%% Merge data
data_out=[SOURCEXY.TIME lat1 lon1];
frontseat_out=[SOURCEXY.TIME navLat_interp navLon_interp];
counter=0;
for ii=1:length(SOURCEXY.TIME)
    if depth_interp(ii)<=.5
%         data_out(ii,1)=navLat_interp(ii);
%         data_out(ii,2)=navLon_interp(ii);
        data_out(ii,:)=NaN;
        frontseat_out(ii,:)=NaN;
        %counter=counter+1;
    end
end
%% Save new LL values
%save(directory + string(start_time) + '_ACOUSTICLL.mat','data_LL','-v7.3')

%% Plot results
% figure
% plot(data_out(:,2),data_out(:,1),navLon_interp,navLat_interp,jetyakLon_interp,jetyakLat_interp)%,lon1,lat1)

% figure
% plot(data.ACOUSTIC_X,data.ACOUSTIC_Y)

%% Generate .kml file
delta_p=zeros(1,length(SOURCEXY.TIME)-1);
for jj=2:length(SOURCEXY.TIME)
    delta_p(jj-1)=sqrt(((SOURCEXY.ACOUSTIC_SOURCE_X(jj)-SOURCEXY.ACOUSTIC_SOURCE_X(jj-1))^2+(SOURCEXY.ACOUSTIC_SOURCE_Y(jj)-SOURCEXY.ACOUSTIC_SOURCE_Y(jj-1))^2))/(SOURCEXY.TIME(jj)-SOURCEXY.TIME(jj-1));
end
trim=false;
data_out_trim=data_out;
frontseat_out_trim=frontseat_out;
status=zeros(1,length(SOURCEXY.TIME)-1);
for kk=2:length(SOURCEXY.TIME)-1
    if isnan(data_out(kk-1,1))
        trim=false;
    else
        if trim==false || delta_p(kk)>2
            data_out_trim(kk,:)=NaN;
            frontseat_out_trim(kk,:)=NaN;
            last=kk-1;
            trim=true;
        elseif trim==true
            delta_p_last=sqrt(((SOURCEXY.ACOUSTIC_SOURCE_X(kk)-SOURCEXY.ACOUSTIC_SOURCE_X(last))^2+(SOURCEXY.ACOUSTIC_SOURCE_Y(kk)-SOURCEXY.ACOUSTIC_SOURCE_Y(last))^2))/(SOURCEXY.TIME(kk)-SOURCEXY.TIME(last));
            if delta_p_last>2
                data_out_trim(kk,:)=NaN;
                frontseat_out_trim(kk,:)=NaN;
            else
                trim=false;
            end
        end
    end
    status(kk)=trim;
end

data_out_trim_trim = data_out_trim(all(~isnan(data_out_trim),2),:); % for nan - rows
frontseat_out_trim_trim = frontseat_out_trim(all(~isnan(frontseat_out_trim),2),:); % for nan - rows
%% Plot figure
figure
plot(data_out_trim_trim(:,3),data_out_trim_trim(:,2),frontseat_out_trim_trim(:,3),frontseat_out_trim_trim(:,2))
%% Continue KML
s=geoshape(data_out_trim_trim(:,2),data_out_trim_trim(:,3));
filename_kml=directory + string(start_time) + "_ACOUSTIC.kml";
kmlwrite(filename_kml,s);

s1=geoshape(jetyakLat_interp,jetyakLon_interp);
filename_kml=directory + string(start_time) + "_JETYAK.kml";
kmlwrite(filename_kml,s1);
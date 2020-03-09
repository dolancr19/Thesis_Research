function [ xUTM, yUTM ] = northriverxy2utm( x, y )
%NORTHRIVERXY2UTM Convert from North River xy curvilinear coordinates to
%UTM coordinates
%   [xUTM,yUTM] = NORTHRIVERXY2UTM(x,y) converts the streamwise coordinates
%   used for the North River to the appropriate UTM coordinates.

%% Load North River coordinates

load('north_river_xy.mat')
dlatb = gradient(latb);
dlonb = gradient(lonb);
thetab = atan2d(dlatb,dlonb);
[xNR,~] = northriver_xy(latb,lonb);

%% Find lat & lon for x

latx = interp1(xNR,latb,x);
lonx = interp1(xNR,lonb,x);
thetax = interp1(xNR,thetab,x);

%% Create map structure

zone = utmzone(nanmean(latb),nanmean(lonb));
ellipsoid= utmgeoid(zone);
mstruct = defaultm('utm');
mstruct.zone = zone;
mstruct.geoid = ellipsoid;
mstruct = defaultm(mstruct);

%% Convert to UTM then adjust by y

[xUTM,yUTM] = mfwdtran(mstruct,latx,lonx);
xUTM = xUTM - y.*sind(thetax);
yUTM = yUTM + y.*cosd(thetax);

end


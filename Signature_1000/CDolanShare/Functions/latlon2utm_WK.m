function [ xutm, yutm ] = latlon2utm_WK( lat, lon)
%LATLON2UTM Projects latitude and longitude coordinates to UTM coordinates
%   This function takes a point of interest to determine the UTM zone and
%   from that zone projects the UTM coordinates. 

p1 = [nanmean(lat), nanmean(lon)];
z1 = utmzone(p1); % utm zone at pod

[ellipsoid, estr] = utmgeoid(z1); % obtain suggested ellipsoid vector and name

utmstruct = defaultm('utm');
utmstruct.zone = z1;
utmstruct.geoid = ellipsoid;
utmstruct = defaultm(utmstruct);

[xutm, yutm] = mfwdtran(utmstruct, lat ,lon);
end

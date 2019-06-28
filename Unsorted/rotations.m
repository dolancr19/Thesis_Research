ele = 30*pi/180;

%azi = -14.4228*(pi/180);
azi = -111*(pi/180);
%pitch = -5.217*(pi/180);
%roll = 3.166*(pi/180);
%hdg = -186.1*(pi/180);

pitch = -10*(pi/180);
roll = 0*(pi/180);
hdg = 186*(pi/180);

rotation_pitch = [1 0 0;0 cos(pitch) sin(pitch);0 -sin(pitch) cos(pitch)];
rotation_roll = [cos(roll) 0 -sin(roll);0 1 0;sin(roll) 0 cos(roll)];
rotation_hdg = [cos(hdg) -sin(hdg) 0;sin(hdg) cos(hdg) 0;0 0 1];

coord = [sin(ele)*cos(azi);sin(ele)*sin(azi);cos(ele)];
coord1 = rotation_pitch*coord;
coord2 = rotation_roll*coord1;
coord3 = rotation_hdg*coord2;

%ele_rot = ((acos(coord3(3)))/sqrt(coord3(1)^2+coord(2)^2+coord(3)^2))*(180/pi)
ele_rot = 180/pi*acos(coord3(3))
azi_rot = (atan(coord3(2)/coord3(1)))*(180/pi)
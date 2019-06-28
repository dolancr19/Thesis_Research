clear variables

spd_current_input = .3;
dir_current_input = 45-90;

stw1 = .5;
spd_error1 = .1;
hdg1 = 000-90;

rotation_hdg = [cosd(hdg1) -sind(hdg1) 0;sind(hdg1) cosd(hdg1) 0;0 0 1];

sog1 = sqrt(((stw1*cosd(hdg1)+spd_current_input*cosd(dir_current_input))^2 + ((stw1*sind(hdg1)+spd_current_input*sind(dir_current_input))^2)));
cog1 = atand((stw1*sind(hdg1)+spd_current_input*sind(dir_current_input))/(stw1*cosd(hdg1)+spd_current_input*cosd(dir_current_input)));

%sog1 = 4;
%cog1 = 000-90;

delta_x1 = (stw1+spd_error1)*cosd(hdg1)+sog1*cosd(cog1);
delta_y1 = (stw1+spd_error1)*sind(hdg1)+sog1*sind(cog1);

stw2 = .5;
spd_error2 = .1;
hdg2 = 180-90;

sog2 = sqrt(((stw2*cosd(hdg2)+spd_current_input*cosd(dir_current_input))^2 + ((stw2*sind(hdg2)+spd_current_input*sind(dir_current_input))^2)));
cog2 = atand((stw2*sind(hdg2)+spd_current_input*sind(dir_current_input))/(stw2*cosd(hdg2)+spd_current_input*cosd(dir_current_input)));

delta_x2 = (stw2+spd_error2)*cosd(hdg2)+sog2*cosd(cog2);
delta_y2 = (stw2+spd_error2)*sind(hdg2)+sog2*sind(cog2);

%stw3 = 5;
%spd_error3 = .5;
%hdg3 = 90-90;

%sog3 = sqrt(((stw2*cosd(hdg2)+spd_current_input*cosd(dir_current_input))^2 + ((stw2*sind(hdg2)+spd_current_input*sind(dir_current_input))^2)));
%cog3 = atand((stw2*sind(hdg2)+spd_current_input*sind(dir_current_input))/(stw2*cosd(hdg2)+spd_current_input*cosd(dir_current_input)));

%delta_x3 = (stw3+spd_error3)*cosd(hdg3)+sog3*cosd(cog3);
%delta_y3 = (stw3+spd_error3)*sind(hdg3)+sog3*sind(cog3);

delta_x = (delta_x1+delta_x2)/2;
delta_y = (delta_y1+delta_y2)/2;

spd_current = sqrt(delta_x^2+delta_y^2)
dir_current = 90+atand(delta_y/delta_x)
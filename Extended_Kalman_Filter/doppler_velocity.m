clear variables
f0=10000; %transmitted frequency in Hz
fd=1; %measured Doppler shift in Hz
c=1500; %speed of sound in m/s
theta=60; %relative angle of source in degrees

%calculate speed through the water using the Doppler shift
%and the relative angle between the receiver and source. 
stw=(fd*c)/(f0*cosd(theta));
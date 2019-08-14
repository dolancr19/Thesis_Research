% script to send usbl pings over the four channel card with built
% in phase shifts
addpath('./jim');
band = 1;
if (band == 1)
  fs = 80000;

% PSK probe Parameters
  fb_probe = 4000;
  fc_probe = 9760;
elseif (band == 2)
  fs = 80000;

% PSK probe Parameters
  fb_probe = 4000;
  fc_probe = 14880;
elseif (band == 3)
  fs = 80000;

% PSK probe Parameters
  fb_probe = 4000;
  fc_probe = 25120;
end
direction = 0;
probe_length = 40/fb_probe ;
silence = 0.100;
nreps = 5;
prepkt = [zeros(floor((silence)*fs),1)];
if direction == 0
  fd1 = fc_probe-fb_probe/2 ;
  fd2 = fc_probe+fb_probe/2 ;
  fb1 = -fb_probe/2;
  fb2 = fb_probe/2;
else
  fd1 = fc_probe+fb_probe/2 ;
  fd2 = fc_probe-fb_probe/2 ;
  fb1 = fb_probe/2;
  fb2 = -fb_probe/2;
end

[y,t]=lfm0(fd1,fd2,probe_length,fs,0);

c = 1430;
lambda = c/fc_probe;
vdist = 0.03;
azimuth0_deg = 45;
elevation0_deg = 0;
azimuth0_rad = (pi/180)*azimuth0_deg;
elevation0_rad = (pi/180)*elevation0_deg;
normal_vector = -[tan(azimuth0_rad) tan(elevation0_rad) 1]';
unit_normal_vector = normal_vector/sqrt(normal_vector'*normal_vector);
v = [0 vdist; vdist vdist; vdist 0 ; 0,0 ]-0.015;
iv = delay_and_phase_shift_array_signal(v, fc_probe, unit_normal_vector, iv1, c_ms, fb_Hz, snr_dB);

dirvec = [ 0.5 0.5 ]';
xunitval = dirvec(1);
yunitval=dirvec(2);
tmp = xunitval*xunitval + yunitval*yunitval ;
if tmp < 1.0
  zunitval = sqrt(1 - tmp);
else
  zunitval = 0.0 ;
end;
az = atan2(yunitval,xunitval) *(180/pi)
el = atan(zunitval/sqrt(xunitval^2 + yunitval^2)) * (180/pi)
az = asin(yunitval/sqrt(xunitval^2 + yunitval^2)) * (180/pi)
numsensors = max(length(v));
T = 1/(fs);
Tcyc = 1/fc_probe;
k = 2*pi/lambda;
 clear phase_shift time_shift y_pb
for i=1:numsensors    
  tau = (v(i,:)*dirvec)*(1/c);
  phase_shift(i) = rem(tau,T) * fc_probe;
  time_shift(i) = rem(tau,T);
  [y,t]=lfm0(fd1,fd2,probe_length,fs);
  y = y.' ;
  y_pb(:,i) = [prepkt;y];
end;
phase_shift*180/pi
len = size(y_pb,1);
y_multi = zeros(len*nreps,numsensors);
for j=1:nreps
  y_multi((j-1)*len+1:j*len,:) = y_pb;
end


out_multi(real(y_multi),fs);
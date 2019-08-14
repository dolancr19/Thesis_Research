%given an azimuth and elevation, calculate the phase difference
%between array elements
% used for testing usbl code
fc = 9760; % carrier
c = 1430;
fs = 80e3;
lambda = c/fc;
count = 1;

%for kk=-1:0.5:1
  %dirvec = [0.1  kk]'
  dirvec = [0 0.5]'
  z = sqrt(1-(dirvec(1).^2+dirvec(2).^2));
  az = atan2(dirvec(1),z)*180/pi;
  el = atan2(dirvec(2),z)*180/pi;
  az_el(count,:) = [az el];

  array_spacing = 0.03; % meters
  v = [0 array_spacing;0 0;  array_spacing 0; array_spacing  array_spacing ];
  numsensors = max(length(v));
  T = 1/(fs);
  Tcyc = 1/fc;
  k = 2*pi/lambda;
  
  for i=1:numsensors
    
    tau(i) = (v(i,:)*dirvec)*(1/c);
    phase_shift(i) = -rem(tau(i),T) * fc;
    time_shift(i) = -rem(tau(i),T);
  end;
  phshift_um2(count,:) = phase_shift*32767;
  count = count+1;
%end

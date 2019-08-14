% function to create array data for one sinusoid and a linear array
% with arbitrary element spacing. Input includes:
%       x               2-d vector of sensor spacings on a plane
%                       [x1 y1; x2 y2; ... xn yn]
%       freq            signal frequency
%       bw              bandwidth (signal time is 1/bw)
%       a	        [ax ay] steering vector to source
%       snr             signal to noise ratio (db)
%       sigtime         length of output (sec)
%       c               sound speed
%       fs              sample freq

%       [y,sig] = arraysig(x,freq,bw,a,snr,sigtime,c,fs)

function [y,sig] =arraysig(x,freq,bw,a,snr,sigtime,c,fs)

  numsensors = max(length(x));
  T = 1/fs;
  Tcyc = 1/freq;
  numpoints = floor(sigtime/T)
  pulse_width = 1/bw;
  noise_mag = 10^(-snr/20);

  lambda = c/freq;
  k = 2*pi/lambda;

  for i=1:numsensors

    tau = (x(i,:)*a)*(1/c);
    delpts = fix(tau/T);
    phase_shift = -rem(tau,T) * freq;
    time_shift = -rem(tau,T);
    sig = sin(2*pi*freq*(time_shift+(0:T:pulse_width)));
    sigsize = round(max(length(sig)));
    sigstart = round(numpoints/2 - sigsize/2);
    fprintf('tau %f, delpts %d, phase_shift %f\n',tau,delpts,phase_shift);

    y(:,i) = noise_mag * randn(numpoints,1);
    y(:,i) = y(:,i) - mean(y(:,i));
    y(sigstart+delpts:sigstart+delpts+sigsize-1,i) =  y(sigstart+delpts:sigstart+delpts+sigsize-1,i) + sig';

  end;






% lfm0.m	L. Freitag	WHOI, July 1992
% function to generate linear fm sweep
% Inputs:	f1	start freq
%		f2	end freq
%		T	time duration
% 		fs	sample frequency
% [y,t] = lfm0(f1,f2,T,fs)

% $Id: lfm0.m 5530 2009-12-11 17:16:54Z jim $

function  [y,t] = lfm0(f1,f2,T,fs,verbose)
 
  finv = 1/fs;
  t=0:finv:(T-finv);
  alpha = 2*pi*(f2-f1)/T;
  wo = 2*pi*f1;
  phi = t.*((0.5*alpha*t)+wo);
  i = sqrt(-1);
  y = exp(i*phi);
 
  if verbose
    Y = fft(y);
    f_axis = (0:length(y)-1)*(fs/length(y));
    plot(f_axis,10*log10(abs(Y)));
   
  end;
 

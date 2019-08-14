% mpingrsp.m
%
% Ping Responder demo script
%
% Sandipa Singh Aug. 2004
% $Id: mpingrsp.m 1784 2005-12-08 16:19:47Z sandipa $
%

  % run this script to issue a $CCRSP cmd to serial port, com1,every 5 sec.
  % then send 25kHz probe to DAQ pad's output channel DAC0 out
% Finally open DAQ pad's input ch. ACH1, to listen to incoming reply.
% Hook up ACH0 to DAC0 out, so as to figure out travel times.
% Entire setup works with hardcoded tat of 50 ms.


port = 'com1';
baudrate = 19200;
Nmax   = 8000;		% max # of samples: 80kHz * 100ms = 8000 samples
nruns = 50; % number of times to run the script


txsig = 1; %25000
rxsig = 1; %25000
timeout = 10000; %ms
tat = 0.05;%50ms

% create fm sweep
fc = 25000;
%fc = 14880;
fs = 80000;
fb = 4000;
nnew=2;
probe_length = 10e-3 ;
t = 1; % seconds to record
exttrig = 0;

[code,fmprobe2]=fm_filt_gen(fb,probe_length,fc,fs,0,'hanning',[],0) ; 
fmpb = b2p(code,[fb,fc]./fs);
fmpb = fmpb./max(fmpb);
y = [fmpb];

% open and configure serial port
sp = serial(port);
set(sp, 'BaudRate', baudrate);
set(sp, 'Terminator', 'CR/LF');
set(sp, 'InputBufferSize', Nmax*2);	% bytes

for i=1:nruns
i
pause(5) % pause 10 sec between issuing commands
fopen(sp);

% issue sweep command
cmd = sprintf('$CCRSP,%d,%d,%d', ...
	rxsig,txsig, timeout);
fprintf(sp, cmd);
disp(sprintf('>>>>> %s', cmd));


% read reply to sweep command
[carsp, carsp_nread] = fscanf(sp, '%s');
if ~carsp_nread
	error('Expected $CARSP');
end
disp(sprintf('<<<<< %s', carsp));


%send sweep via daq pad
[tt(i),yrcv] = out_and_in_on_trig(y,fs,tat,1,t,exttrig, fmprobe2);

% close serial port
fclose(sp);

end % for i=1:nruns




% to change Micromodem baud rate:   $CCCFG,BR2,N  where N is:
%
%	0	  2400
%	1	  4800
%	2	  9600
%	3	 19200
%	4	 38400
%	5	 57600
%	6	115200
%	7	230400	(untested with RS-232 isolation circuit!)


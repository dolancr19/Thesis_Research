% test umodem's rsp mode. 
%send out a ping and wait for rcvd. ping
%print out one way travel time
% S. Singh

% create fm sweep
fc = 25000;
%fc = 14880;
fs = 80000;
fb = 4000;
nnew=2;
probe_length = 10e-3 ;

[code,fmprobe2]=fm_filt_gen(fb,probe_length,fc,fs,0,'hanning',[],0) ; 
fmpb = b2p(code,[fb,fc]./fs);
fmpb = fmpb./max(fmpb);
y = [fmpb];
t = 1; % seconds to record
exttrig = 0;

[tt,yrcv] = out_and_in_on_trig(y,fs,0.2,1,t,exttrig, fmprobe2);
if 0
tat=0.05;
% mfd the data
%bb = p2b(yrcv, [fb*nnew,fc]/fs);
[irefm, psfm, dv,p] = mfd1(bb(:,1:2:end), conj(fmprobe2),600);
index = find(psfm(:,1)==-1);
psfm(index,:) = [];
bborig = bb;

bb(2*(index-1)+1:2*index,:) = [];
clear tt2;
for (i=1:size(psfm,1))
     [irefm, psfm2, dv,p] = mfd1(bb(psfm(i,1):end,2*i), conj(fmprobe2),600);
     tt2(i) = psfm2(1,1)/(fb*nnew);
end
tt = (tt2-tat)/2
end
if 0
     time(index) = [];
tb = time/10-psfm(:,1).'
end

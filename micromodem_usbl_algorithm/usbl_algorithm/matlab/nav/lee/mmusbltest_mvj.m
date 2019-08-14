% mmusbltest.m     script for testing beamforming with pulses

% 2018/09/26 20:45:15  MVJ  This version is for a narrow-band USBL.  See usbl_test_jwp_mvj.m



fc = 15000; % carrier
fs = 80000 ; % passband  - as in sample rate.  There must be an analog anti-aliasing filter ahead of this.
bw = 500; % Not a sweep, just defines the width of a const. freq. pulse at fc.
nnew = 2 ;
snr = 20;
c = 1500;
total_time = 0.005;
detv = [1 1] ;
dthresh = 10 ;
pown = 20 ;

lambda = c/fc;

dir_vec = [ 0.5 0.5  ]'

v = [0 0; 0 lambda/2; lambda/2 lambda/2 ; lambda/2 0 ];

% this creates the input signal, appropriately delayed for each element.  But the signal is
% a pulse at constant frequency (freq) for a duration (1/bw), not a sweep.
[y,sig_replica] = twodasig(v, fc, bw, dir_vec, snr, total_time, c, fs);

% passband to baseband deconvolution.  This should yield a pulse absent any noise.
iv = p2b(y,[bw*nnew fc]/fs) ;

% matched filter detection.  detv = [1 1] is the code, i.e. the filter coefficients for a pulse, I guess?
[ire,ps,dv,p] = mfd1(iv,detv,dthresh,pown) ;
[mval,mloc] = max(ire(:,1)) ;

X = ire(mloc,:) ;
R = X'*X ;


[S,xs,ys,A] = td2bf(R,v,fc,c,[-1 1 -1 1],0.05);
[u,az_el] = bf_fast(R,v,fc,c,[-1 1 -1 1],0.05,10);
[maxr,maxrloc] = max(S);
[maxc,maxcloc] = max(maxr);
xloc = maxrloc(maxcloc);
yloc = maxcloc;
solution_vec = [xs(xloc) ys(yloc)]

figure(2) ; 
pcolor(xs,ys,S) ;grid
z = sqrt(1-(dir_vec(1).^2+dir_vec(2).^2));
az = atan(dir_vec(1)/z)*180/pi
el = atan(dir_vec(2)/z)*180/pi

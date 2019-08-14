% mmusbltest.m     script for testing beamforming with pulses


fc = 15000; % carrier
fs = 80000 ; % passband
bw = 5000; 
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

[y,sig_replica] = twodasig(v, fc, bw, dir_vec, snr, total_time, c, fs);

iv = p2b(y,[bw*nnew fc]/fs) ;

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

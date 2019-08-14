% mmusbltest.m     script for testing beamforming with pulses


fc = 9760; % carrier
fs = 80000 ; % passband
bw = 4000; 
nnew = 2 ;
snr = 20;
c = 1430;
total_time = 0.005;
detv = [1 1] ;
dthresh = 10 ;
pown = 20 ;

lambda = c/fc;
vdist = 0.03;
dir_vec = [ 0 1  ]'

v = [0 vdist; 0 0; vdist 0 ; vdist,vdist ]-0.015;

[y,sig_replica] = twodasig(v, fc, bw, dir_vec, snr, total_time, c, fs);

iv = p2b(y,[bw*nnew fc]/fs) ;

[ire,ps,dv,p] = mfd1(iv,detv,dthresh,pown) ;
[mval,mloc] = max(ire(:,1)) ;

X = ire(mloc,:) ;
diff(phase(X))*180/pi
R = X'*X ;
% to add a known phase difference to um2 ire calcs, use this
cumsum(diff(phase(X))*32767/(2*pi))

[S,xs,ys,A] = td2bf(R,v,fc,c,[-1 1 -1 1],0.05);
[u,az_el] = bf_fast(R,v,fc,c,[-1 1 -1 1],0.05,10)
[maxr,maxrloc] = max(S);
[maxc,maxcloc] = max(maxr);
xloc = maxrloc(maxcloc);
yloc = maxcloc;
solution_vec = [xs(xloc) ys(yloc)]
xunitval = solution_vec(1);
yunitval=solution_vec(2);
tmp = xunitval*xunitval + yunitval*yunitval ;
if tmp < 1.0
  zunitval = sqrt(1 - tmp);
else
  zunitval = 0.0 ;
end;
if xunitval == 0.0
  xunitval = eps;
end;
if yunitval == 0.0
  yunitval == eps;
end;

if xunitval >= 0
  if yunitval >= 0
    az = atan(abs(xunitval)/abs(yunitval))*(180/pi);
  end;
  if yunitval < 0
    az = 90 + atan(abs(yunitval)/abs(xunitval))*(180/pi);
  end;
end;
if xunitval < 0
  if yunitval >= 0
    az = 270 + atan(abs(yunitval)/abs(xunitval))*(180/pi);
  end;
  if yunitval < 0
    az = 180 + atan(abs(xunitval)/abs(yunitval))*(180/pi);
  end;
end;

if real(zunitval) == 0.0
  zunitval = 0.0;
end;

el = atan(zunitval/sqrt(xunitval^2 + yunitval^2)) * (180/pi);
az,el


figure(2) ; 
pcolor(xs,ys,S) ;grid

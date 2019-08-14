DIAM = 0.300; % m  Torqeedo 1901
RHO = 1030;
CNTS2RPM = 1/(12*3)*60;  % this from dualddt tests apparently confirmed with dyno.

Ktau = 0.601 ;% (N m/A)                Torque constant and voltage constant are numerically identical.
Kb = 0.601; % (V s/rad)      Back emf constant
Rm = 0.794; % (Ohm)          Terminal resistance


% Estimate Zww.  This is from clio013, a mid-depth ascent between samples and the descent.
% a more involved analysis is in nui015/eng/ascent_descent.m  This is derived from that.
TWINA =  1.0e+09 *[   1.540358882387986   1.540358969129355];
TWIND =  1.0e+09 *[   1.540265743288343   1.540266444511677];
TWINZ = 1.0e+09 *[   1.540181836127314   1.540185717770469];


cssd = sextt(css,TWIND);
cssa = sextt(css,TWINA);
cssz = sextt(css,TWINZ);

% NUI fit underpredicts current during descent.
% inclined to believe this over Sentry fit because wake fraction should reduce thruster effectiveness on vehicle
% over open water flume results.  That's for AV.  Prediction is pretty spot on for PV/SV in Q1.
pKT = [-0.1075   -0.2487    0.3742];  % from flume data for NUI torqeedo.
pKQtau = [0.0320   -0.0337    0.0611];

% zero-order reverse-bollard fits.  Only NUI thruster data available here.
pKTr = [-0.1820];
pKQtaur = [-0.0509];

elm.PV = elm0; % assignment is arbitrary.  I don't know which is which.
elm.SV = elm1; 

fn={'pv','sv'};
for n=1:length(fn)

  s = elm.(upper(fn{n}));
  s.w = interp1(css.t,css.vel(:,3),s.t);
  s.rps = -s.en_vel*CNTS2RPM/60;  % change sign so that positive is pushing down. 
  s.J = s.w./s.rps/DIAM;
  s.irev = s.rps < 0;
  s.KT = polyval(pKT,s.J);
  s.KT(s.irev) = polyval(pKTr,s.J(s.irev));
  s.T = s.KT.*RHO.*s.rps.^2.*DIAM^4;
  s.KQtau = polyval(pKQtau,s.J);
  s.KQtau(s.irev) = polyval(pKQtaur,s.J(s.irev));
  s.Im = -s.KQtau.*RHO.*s.rps.^2.*DIAM^5;  % use this to cross-check against measured. 
  s.Q = s.Im*Ktau;
  s.Pprop = s.T.*s.w;
  s.Pshaft = abs(s.Q.*s.rps*2*pi);
  s.Vemf = s.rps*2*pi.*Kb;  % does not include voltage from winding resistance.
  s.Pm = abs(s.Vemf.*s.Im) + s.Im.^2*Rm;
  s.Psupply = s.Pm;
  
  % all of the above is predicted.  Now do measured.
  s.Im_obs = s.current;
  s.Q_obs = s.Im_obs*Ktau;
  s.Pshaft_obs = abs(s.Q_obs.*s.rps*2*pi);
  s.Pm_obs = abs(s.Vemf.*s.Im_obs) + s.Im_obs.^2*Rm;
  s.Psupply_obs = s.Pm_obs;
  
  eval(sprintf('%s = s;',fn{n}));
  eval(sprintf('%sd = sextt(%s,TWIND);',fn{n},fn{n}));
  eval(sprintf('%sa = sextt(%s,TWINA);',fn{n},fn{n}));
  eval(sprintf('%sz = sextt(%s,TWINZ);',fn{n},fn{n}));
end

% cross-check predicted and actual current.
% looks like prediction is pretty good during descent, pretty bad, but conservative, during ascent.
% This makes sense because we don't have reliable reverse data at speed, only bollard.
figure(2); clf reset;
subplot(211)
plot(pv.t,pv.Im,pv.t,pv.current);
legend('Predicted','Observed');
ylabel('PV Current (A)');
ylim([-30 30]);
subplot(212)
plot(sv.t,sv.Im,sv.t,sv.current);
ylabel('SV Current (A)');
ylim([-30 30]);

% vehicle drag
% this plot says 350 N down, -350 up, and about 100 N holding depth.  The up data is wrong because the vehicle is
% moving fast.  The thrusters are definitely generating a fair bit less thrust than they would at bollard.
figure(5); clf reset;
pv = sunique(pv);
sv = sunique(sv);
v.t = pv.t;
v.Tpv = interp1(pv.t,pv.T,v.t);
v.Tsv = interp1(sv.t,sv.T,v.t);
v.Z = v.Tpv+v.Tsv;
plot(v.t,v.Z);
ylabel('Z (N)');

% Because of uncertainty in ascent data, will only use descent and hold data to compute Zww.
% now compute fits to nominally steady-state ascent and descent conditions.
wd = robustfit(cssd.t-cssd.t(1),-cssd.pos(:,3));
wd = wd(2);

Zd = robustfit(pvd.t-pvd.t(1),pvd.T) + robustfit(svd.t-svd.t(1),svd.T);
Zd = Zd(1);


Zz = robustfit(pvz.t-pvz.t(1),pvz.T) + robustfit(svz.t-svz.t(1),svz.T);
Zz = Zz(1);

Zww = (Zd-Zz)/wd^2;
fprintf(1,'Zww computed from descent plus depth hold data (for ballast) is: %.1f N s^2/m^2\n',Zww);
fprintf(1,'Caveat!  Vehicle was in unstable speed regime during descent!\n');

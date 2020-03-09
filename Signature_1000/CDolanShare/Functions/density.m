function [ rho ] = density( s, t )
%DENSITY Returns the density of seawater given salinity and temperature
%   rho = DENSITY(s,t) calculates the density of seawater in kg/m3 given
%   salinity (s) in PSU and temperature (t) in C. The equation of state is
%   calculated using the empirical formulation given in Adrian E. Gill's
%   Atmosphere-Ocean Dynamics (1983).
%
%   A. P. Garcia

rhow = 999.842594 + 6.793952E-2*t - 9.095290E-3*t.^2 + 1.001685E-4*t.^3 ...
    - 1.120083E-6*t.^4 + 6.536332E-9*t.^5;

rho = rhow + s.*(0.824493 - 4.0899E-3*t + 7.6438E-5*t.^2 ...
    - 8.2467E-7*t.^3 + 5.3875E-9*t.^4) ...
    + s.^(3/2).*(-5.72466E-3 + 1.0227E-4*t - 1.6546E-6*t.^2) ...
    + 4.8314E-4*s.^2;

end
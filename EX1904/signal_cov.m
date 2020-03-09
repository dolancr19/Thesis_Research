function [Sx]=signal_cov(x,n_snap,n_sensor)
%% Prepare workspace
% clear variables
% clc

%% Create signal covariance matrix
Sx=zeros(n_sensor);

for ii=1:n_snap
    cov=x(:,ii)*x(:,ii)';
    Sx=Sx+cov;
end
Sx=Sx/n_snap;
end
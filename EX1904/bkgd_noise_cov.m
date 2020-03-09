function [Sn]=bkgd_noise_cov(n,n_snap,n_sensor)
%% Prepare workspace
% clear variables
% clc

%% Create background noise covariance matrix
Sn=zeros(n_sensor);

for ii=1:n_snap
    cov=n(:,ii)*n(:,ii)';
    Sn=Sn+cov;
end
Sn=Sn/n_snap;
end
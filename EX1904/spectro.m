clear variables
clc

%% Common variables
%load combined1.mat
load D:\Documents\Thesis_Research\EX1904\Processed_Data\data55.mat
fs=45000;
steps=length(data);
begin=364500;
last=378000;

n_snap=last-begin; 
n_sensor=2;
f0=10000; %%% the frequency to beamform
c=1500; %%% water sound speed
%d=c/f0/2; %%% array spacing (?)

%p_=(0:n_sensor-1).'*d;
p_=[.019;-.019];
%% Plot frequency domain

% for ii=2:5
%     y=data(:,ii);
%     y_bp=bandpass(y,[7500 12500],45000);
%     y_f=fft(y_bp,steps);
%     f=(0:steps-1)*fs/steps;
% 
%     figure
%     plot(f,10*log10(abs(y_f).^2))
%     xlabel('Frequency (Hz)')
%     ylabel('Amplitude')
%     xlim([0 fs/2])
% end

%% Plot spectrogram
win=rectwin(512);
noverlap=length(win)/2;
nfft=length(win);
%y0=zeros(round(steps/256),n_sensor);
for ii=2:5
    y=data(:,ii);
    y_bp=bandpass(y,[7500 12500],45000);
    [s,f,t]=spectrogram(y,win,noverlap,nfft,fs);
    if ii==2
        [~,ind_f]=min(abs(f-f0));
    end
    y0(:,ii)=s(ind_f,:);
%     figure
%     imagesc(t,f,10*log10(abs(s).^2))
%     ylim([0 fs/2])
%     axis xy
%     xlabel('Time (s)')
%     ylabel('Frequency (Hz)')
%     colormap('winter')
end

%load 'D:\Documents\Thesis_Research\EX1904\background_noise_7_31.mat'
%% Simu param

% SNR=0;   %%% in dB
% 
% %%% Signal of interest
% theta_t=95;  %%% target position (degrees) 
% Ft=5; %%% target amplitude [we assume this is unknown]
% 
% %%% Interfering signal 
% theta_interf=125; %%% interference position 
% Fi=50; %%% interference amplitude 



%% Create noise
% mu = zeros(1,n_sensor);
% sigma = eye(n_sensor);
% 
% noise_r = mvnrnd(mu,sigma,n_snap);
% noise_i = mvnrnd(mu,sigma,n_snap);
% noise=noise_r+1i*noise_i;
% 
% noise=noise.';



%% Create target (source) and interference signals
% v_t=exp(1i*cosd(theta_t)*p_*2*pi*f0/c);
% v_interf=exp(1i*cosd(theta_interf)*p_*2*pi*f0/c);
% 
% v_t_rep=Ft*repmat(v_t,1,n_snap);
% v_interf_rep=Fi*repmat(v_interf,1,n_snap);
% 
% %%%% Add random phase shift, because we don't know where we start the
% %%%% snapshot
% for ii=1:n_snap
%     a= 2*pi.*rand(1,1)-pi;  
%     b= 2*pi.*rand(1,1)-pi;  
%     v_t_rep(:,ii)=v_t_rep(:,ii)*exp(1i*a);
%     v_interf_rep(:,ii)=v_interf_rep(:,ii)*exp(1i*b);
% end



%% Adjust SNR of ambient noise
% Ps=sum(abs(v_t).^2);
% Pn_obj=Ps*10^(-SNR/10);
% 
% Pn=mean(sum(abs(noise).^2,1));
% noise_ok=noise/sqrt(Pn)*sqrt(Pn_obj);
% 
% %%%% verif
% Pn_ok=mean(sum(abs(noise_ok).^2,1));
% SNR_xp=10*log10(Ps/Pn_ok);

%% Create noisy signal
% x_in=v_t_rep+v_interf_rep+noise_ok;

%% CBF Weights
w_cbf=ones(n_sensor,1)/n_sensor;
%% MVDR weights
x_in=data.';
n=[x_in(2,fs*7.5:fs*7.6);x_in(4,fs*7.5:fs*7.6)];
[Sn]=bkgd_noise_cov(n,length(n),n_sensor);
% middle=round(length(n)/2);
% n_t=n(:,middle);


%% MPDR weights
%x=combined(1:n_snap,:)';

x_samp=[x_in(2,begin:last);x_in(4,begin:last)];
x=[data(371250,2);data(371250,4)];
[Sx]=signal_cov(x_samp,n_snap,n_sensor);



%% Beampatterns

angle=0:.1:180;

b_cbf=zeros(1,length(angle));
b_mvdr=zeros(1,length(angle));
b_mpdr=zeros(1,length(angle));


for ii=1:length(angle)
    k_s=(2*pi*f0/c)*cosd(angle(ii));
    v_s=exp(1i*k_s*p_);
    
    w_mvdr=inv(Sn)*v_s/(v_s'*inv(Sn)*v_s);
    w_mpdr=inv(Sx)*v_s/(v_s'*inv(Sx)*v_s);
    
    b_cbf(ii)=w_cbf'*x;
    b_mvdr(ii)=w_mvdr'*x;
    b_mpdr(ii)=w_mpdr'*x;
    
end

figure
plot(angle,abs(b_cbf))
 hold on
 plot(angle,(abs(b_mvdr)))
 hold on
 plot(angle,(abs(b_mpdr)))
 legend('Conventional', 'MVDR', 'MPDR')

% %% Compute and plot Y
% 
% Y_cbf=zeros(1,n_snap);
%  Y_mvdr=zeros(1,n_snap);
%  Y_mpdr=zeros(1,n_snap);
% 
% for ii=1:n_snap
%     y_cbf=w_cbf'*x_samp(:,ii);
%      y_mvdr=w_mvdr'*x_samp(:,ii);
%      y_mpdr=w_mpdr'*x_samp(:,ii);
%     
%     Y_cbf(ii)=y_cbf;
%      Y_mvdr(ii)=y_mvdr;
%      Y_mpdr(ii)=y_mpdr;
% end
% 
% % figure
% % plot(1:n_snap,abs(Y_cbf))
% %  hold on
% %  plot(1:n_snap,abs(Y_mvdr))
% %  hold on
% %  plot(1:n_snap,abs(Y_mpdr))
% %  legend('Conventional', 'MVDR', 'MPDR')
% figure
% imagesc(angle,t,abs(Y_cbf))
% figure
% imagesc(angle,t,abs(Y_mvdr))
% figure
% imagesc(angle,t,abs(Y_mpdr))
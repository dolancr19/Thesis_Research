%% Prepare workspace
clear variables;
close all;
clc;

%% Define contants
iterations=10;
freq=1; %Number of cycles per second
steps=3600*freq; %# of filter steps desired
track0=build_track(freq,steps); %build simulated values for speed through the water and heading
speeds=[.1];
headings=[0];

v_k_full=randn(4,steps);
set0=.25;
drift0=90-60;
variation=0; % zero for constant current, non-zero defines band for oscillating current
% current=build_current(steps, set0, drift0, variation);

% data_full=zeros(27,steps,iterations); % NCV_C
% mu_full=zeros(6,steps,iterations); % NCV_C
% mu_bar=zeros(6,steps); % NCV_C
data_full=zeros(29,steps,iterations); % NCV_C_bias
mu_full=zeros(7,steps,iterations); % NCV_C_bias
mu_bar=zeros(7,steps); % NCV_C_bias
% data_full=zeros(29,steps,iterations); % NCV_bias
% mu_full=zeros(7,steps,iterations); % NCV_bias
% mu_bar=zeros(7,steps); % NCV_bias
epsilon_full=zeros(iterations,steps);
epsilon_bar=zeros(1,steps);
epsilon_v_full=zeros(iterations,steps);
epsilon_v_bar=zeros(1,steps);





%Define boundary conditions
x0_pos=10;
y0_pos=20;
stw0=.50;
hdg0=90-30;
%% Execute filter
% Add bias for testing
for oo=1:length(speeds)
    track.stw=track0.stw+speeds(oo);
    stw_bias=speeds(oo);
    for pp=1:length(headings)
        track.hdg=track0.hdg+headings(pp);

for ii=1:iterations
    % Uncomment below for different current for each run
%     set_lower=0;
%     set_upper=.5;
%     drift_lower=-180;
%     drift_upper=180;
% 
%     set0=(set_upper-set_lower)*rand(1,1)+set_lower;
%     drift0=(drift_upper-drift_lower)*rand(1,1)+drift_lower;
    current=build_current(steps, set0, drift0, variation);

%     [data,epsilon,epsilon_v,mu,K_out, H_out,P_plus_out,P_minus_out]=NCV(freq, steps, x0_pos, y0_pos, stw0, hdg0,track, current);
%     [data,epsilon,epsilon_v,mu,K_out, H_out,P_plus_out,P_minus_out]=NCV_bias(freq, steps, x0_pos, y0_pos, stw0, hdg0,track, current,track0);
%     [data,epsilon,epsilon_v,mu,K_out, H_out,P_plus_out,P_minus_out]=NCV_C(freq, steps, x0_pos, y0_pos, stw0, hdg0,track, current,v_k_full,track0);
%     [data,epsilon,epsilon_v,mu,K_out, H_out,P_plus_out,P_minus_out]=NCV_C_bias(freq, steps, x0_pos, y0_pos, stw0, hdg0,track, current,v_k_full,track0,stw_bias);
    [data,epsilon,epsilon_v,mu,K_out, H_out,P_plus_out,P_minus_out]=NCV_C_bias_hybrid(freq, steps, x0_pos, y0_pos, stw0, hdg0,track, current,v_k_full,track0,stw_bias);
    data_full(:,:,ii)=data;
    epsilon_full(ii,:)=epsilon;
    epsilon_bar=epsilon_bar+epsilon;
    epsilon_v_full(ii,:)=epsilon_v;
    epsilon_v_bar=epsilon_v_bar+epsilon_v;
    mu_full(:,:,ii)=mu;
    mu_bar=mu_bar+mu;
end

epsilon_bar=epsilon_bar/iterations;
epsilon_v_bar=epsilon_v_bar/iterations;
mu_bar=mu_bar/iterations;
save(string(speeds(oo))+'_'+string(headings(pp))+'.mat','data_full','epsilon_bar','epsilon_v_bar','mu_bar','-v7.3')
    end
end

%% Plot results


chi_sq_upper1=chi2inv(0.995,6*iterations)/iterations;
chi_sq_upper2=chi2inv(0.995,4*iterations)/iterations;
chi_sq_lower1=chi2inv(0.005,6*iterations)/iterations;
chi_sq_lower2=chi2inv(0.005,4*iterations)/iterations;

figure(1)
plot(epsilon_bar)
yline(chi_sq_upper1,'--b');
yline(chi_sq_lower1,'--b');
%title('10 Run Normalized Estimation Error Squared (NEES)')
xlabel('Filter Step')
xlim([0 steps])
ylim([0 20])
ylabel('NEES')
set(gca,'fontsize',16);

figure(2)
plot(epsilon_v_bar)
yline(chi_sq_upper2,'--b');
yline(chi_sq_lower2,'--b');
%title('10 Run Normalized Innovation Squared (NIS)')
xlabel('Filter Step')
xlim([0 steps])
ylabel('NIS')
set(gca,'fontsize',16);

figure(3)
% plot(1:steps,data_full(17,:,1)) % NCV_C
% plot(1:steps,data_full(21,:,1)) % NCV_C_bias
plot(1:steps,data_full(19,:,1)) % NCV_bias
%title('Error Covariance Matrix Norm vs Filter Step')
xlabel('Filter Step')
xlim([0 steps])
ylabel('Error Covariance Matrix Norm')
set(gca,'fontsize',16);

figure(4)
reset(gca)
yyaxis left
% plot(1:steps,data_full(18,:,1)) % NCV_C
% plot(1:steps,data_full(22,:,1)) %NCV_C_bias
plot(1:steps,data_full(20,:,1)) % NCV_bias
%title('Filter Error vs Filter Step')
xlabel('Filter Step')
xlim([0 steps])
ylabel('Filter Error (m)')
yyaxis right
% plot(1:steps,data_full(19,:,1)) % NCV_C
% plot(1:steps,data_full(23,:,1)) % NCV_C_bias
plot(1:steps,data_full(21,:,1)) % NCV_bias
ylabel('Filter Error (degrees)')
legend('Position Error','Azimuth Error')
set(gca,'fontsize',16);

figure(5)
reset(gca)
% plot(data_full(7,:,1),data_full(9,:,1),data_full(15,:,1),data_full(16,:,1)) % NCV_C
% plot(data_full(9,:,1),data_full(11,:,1),data_full(19,:,1),data_full(20,:,1)) % NCV_C_bias
plot(data_full(8,:,1),data_full(10,:,1),data_full(17,:,1),data_full(18,:,1)) % NCV_bias

axis equal
%title('Acutal Position vs Filter Output Position')
xlabel('Easting position (m)')
ylabel('Northing position (m)')
legend('Filter Output','Actual Position','location','southeast')
set(gca,'fontsize',16);

figure(6)
reset(gca)
ax1 = subplot(6,1,1); % NCV_C
% ax1 = subplot(8,1,1); %NCV_C_bias
x = 1:steps;
y1 = mu_bar(1,:); 
plot(ax1,x,y1)
axis([0 steps -3 3])
title('Easting Position')
xlabel('Filter Step')
ylabel('NMEE')
set(gca,'fontsize',12);

ax2 = subplot(6,1,2); % NCV_C
% ax2 = subplot(8,1,2); % NCV_C_bias
x = 1:steps;
y2 = mu_bar(2,:); 
plot(ax2,x,y2)
axis([0 steps -3 3])
title('Easting STW')
xlabel('Filter Step')
ylabel('NMEE')
set(gca,'fontsize',12);

ax3 = subplot(6,1,3); % NCV_C
% ax3 = subplot(8,1,3); % NCV_C_bias
x = 1:steps;
y3 = mu_bar(3,:); 
plot(ax3,x,y3)
axis([0 steps -3 3])
title('Northing Position')
xlabel('Filter Step')
ylabel('NMEE')
set(gca,'fontsize',12);

ax4 = subplot(6,1,4); % NCV_C
% ax4 = subplot(8,1,4); % NCV_C_bias
x = 1:steps;
y4 = mu_bar(4,:); 
plot(ax4,x,y4)
axis([0 steps -3 3])
title('Northing STW')
xlabel('Filter Step')
ylabel('NMEE')
set(gca,'fontsize',12);

ax5 = subplot(6,1,5); % NCV_C
% ax5 = subplot(8,1,5); % NCV_C_bias
x = 1:steps;
y5 = mu_bar(5,:); 
plot(ax5,x,y5)
axis([0 steps -3 3])
title('Easting Current')
xlabel('Filter Step')
ylabel('NMEE')
set(gca,'fontsize',12);

ax6 = subplot(6,1,6); % NCV_C
% ax6 = subplot(8,1,6); % NCV_C_bias
x = 1:steps;
y6 = mu_bar(6,:); 
plot(ax6,x,y6)
axis([0 steps -3 3])
title('Northing Current')
xlabel('Filter Step')
ylabel('NMEE')
set(gca,'fontsize',12);

% ax7 = subplot(8,1,7); % NCV_C_bias
% x = 1:steps;
% y7 = mu_bar(7,:); 
% plot(ax7,x,y7)
% axis([0 steps -3 3])
% title('Easting Bias')
% xlabel('Filter Step')
% ylabel('NMEE')
% set(gca,'fontsize',12);
% 
% ax8 = subplot(8,1,8); % NCV_C_bias
% x = 1:steps;
% y8 = mu_bar(8,:); 
% plot(ax8,x,y8)
% axis([0 steps -3 3])
% title('Northing Bias')
% xlabel('Filter Step')
% ylabel('NMEE')
% set(gca,'fontsize',12);

yline(ax1,2.57/sqrt(iterations),'--b');
yline(ax2,2.57/sqrt(iterations),'--b');
yline(ax3,2.57/sqrt(iterations),'--b');
yline(ax4,2.57/sqrt(iterations),'--b');
yline(ax5,2.57/sqrt(iterations),'--b');
yline(ax6,2.57/sqrt(iterations),'--b');
% yline(ax7,2.57/sqrt(iterations),'--b');
% yline(ax8,2.57/sqrt(iterations),'--b');
yline(ax1,-2.57/sqrt(iterations),'--b');
yline(ax2,-2.57/sqrt(iterations),'--b');
yline(ax3,-2.57/sqrt(iterations),'--b');
yline(ax4,-2.57/sqrt(iterations),'--b');
yline(ax5,-2.57/sqrt(iterations),'--b');
yline(ax6,-2.57/sqrt(iterations),'--b');
% yline(ax7,-2.57/sqrt(iterations),'--b');
% yline(ax8,-2.57/sqrt(iterations),'--b');
% %sgtitle('10 Run Normalized Mean Estimation Error (NMEE)')

% NCV_C
% mag=hypot(data_full(11,:,1),data_full(12,:,1)); 
% dir=atan2d(data_full(12,:,1),data_full(11,:,1));
% NCV_C_bias
% mag=hypot(data_full(13,:,1),data_full(14,:,1));
% dir=atan2d(data_full(14,:,1),data_full(13,:,1));
% NCV_C
mag=hypot(data_full(12,:,1),data_full(13,:,1)); 
dir=atan2d(data_full(13,:,1),data_full(12,:,1));

figure(7)
reset(gca)
yyaxis left
plot(1:steps,current.set,1:steps,mag)
%title('Calculated Current vs Filter Step')
xlabel('Filter Step')
ylabel('Magnitude (m/s)')
axis([1 steps 0 2])
yyaxis right
plot(1:steps,90.-current.drift,1:steps,90.-dir)
ylabel('Direction (degrees)')
axis([1 steps -180 180])
set(gca,'fontsize',16);

% bias_mag=hypot(data_full(32,:,1),data_full(33,:,1));
% bias_calc=hypot(data_full(15,:,1),data_full(16,:,1));
% bias_dir=atan2d(data_full(16,:,1),data_full(15,:,1));

% 
% figure(8)
% reset(gca)
% plot(1:steps,bias_mag,1:steps,bias_calc)
% 
% figure(9)
% reset(gca)
% plot(1:steps,track.hdg,1:steps,bias_dir)

figure(10)
reset(gca)
plot(1:steps,data_full(14,:,1))
hold on
yline((stw0+stw_bias)/stw0)
xlim([0 steps])
xlabel('Filter Step')
ylabel('Bias Ratio')
legend('Calculated','Actual')
set(gca,'fontsize',12);

mean(epsilon_bar)
mean(epsilon_v_bar)
% mean(bias_calc)
%% Prepare workspace
clear variables;
clc

%% Define contants
iterations=10;
freq=10; %Number of cycles per second
steps=4000*freq; %# of filter steps desired
track=build_track(freq,steps); %build simulated values for speed through the water and heading

set0=.25;
drift0=60;
variation=0; % zero for constant current, non-zero defines band for oscillating current
% current=build_current(steps, set0, drift0, variation);

data_full=zeros(27,steps,iterations); %matrix for post-processing data
epsilon_full=zeros(iterations,steps);
epsilon_bar=zeros(1,steps);
epsilon_v_full=zeros(iterations,steps);
epsilon_v_bar=zeros(1,steps);
mu_full=zeros(6,steps,iterations);
mu_bar=zeros(6,steps);

%Define boundary conditions
x0_pos=10;
y0_pos=20;
stw0=.50;
hdg0=30;

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
    [data,epsilon,epsilon_v,mu,K_out, H_out,P_plus_out,P_minus_out]=NCV_C(freq, steps, x0_pos, y0_pos, stw0, hdg0,track, current);
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


%% Plot results
chi_sq_upper1=chi2inv(0.995,6*iterations)/iterations;
chi_sq_upper2=chi2inv(0.995,4*iterations)/iterations;
chi_sq_lower1=chi2inv(0.005,6*iterations)/iterations;
chi_sq_lower2=chi2inv(0.005,4*iterations)/iterations;

figure(1)
plot(epsilon_bar)
yline(chi_sq_upper1,'--b');
yline(chi_sq_lower1,'--b');

figure(2)
plot(epsilon_v_bar)
yline(chi_sq_upper2,'--b');
yline(chi_sq_lower2,'--b');

figure
plot(1:steps,data_full(17,:,1))
title('Error Covariance Matrix Norm vs Filter Step')
xlabel('Filter Step')
ylabel('Error Covariance Matrix Norm')

figure
plot(1:steps,data_full(18,:,1))
%axis([0 steps -5 5])
title('Filter Error vs Filter Step')
xlabel('Filter Step')
ylabel('Filter Error (m)')

figure
plot(1:steps,data_full(19,:,1))
%axis([0 steps -5 5])
title('Filter Error vs Filter Step')
xlabel('Filter Step')
ylabel('Filter Error (degrees)')

figure
plot(data_full(7,:,1),data_full(9,:,1),data_full(15,:,1),data_full(16,:,1))
axis equal
title('Acutal Position vs Filter Output Position')
xlabel('x position (m)')
ylabel('y position (m)')
legend('Filter Output','Actual Position','location','southeast')

figure
ax1 = subplot(6,1,1);
x = 1:steps;
y1 = mu_bar(1,:); 
plot(ax1,x,y1)
axis([0 steps -3 3])

ax2 = subplot(6,1,2);
x = 1:steps;
y2 = mu_bar(2,:); 
plot(ax2,x,y2)
axis([0 steps -3 3])

ax3 = subplot(6,1,3);
x = 1:steps;
y3 = mu_bar(3,:); 
plot(ax3,x,y3)
axis([0 steps -3 3])

ax4 = subplot(6,1,4);
x = 1:steps;
y4 = mu_bar(4,:); 
plot(ax4,x,y4)
axis([0 steps -3 3])

ax5 = subplot(6,1,5);
x = 1:steps;
y5 = mu_bar(5,:); 
plot(ax5,x,y5)
axis([0 steps -3 3])

ax6 = subplot(6,1,6);
x = 1:steps;
y6 = mu_bar(6,:); 
plot(ax6,x,y6)
axis([0 steps -3 3])

yline(ax1,2.57/sqrt(iterations),'--b');
yline(ax2,2.57/sqrt(iterations),'--b');
yline(ax3,2.57/sqrt(iterations),'--b');
yline(ax4,2.57/sqrt(iterations),'--b');
yline(ax5,2.57/sqrt(iterations),'--b');
yline(ax6,2.57/sqrt(iterations),'--b');
yline(ax1,-2.57/sqrt(iterations),'--b');
yline(ax2,-2.57/sqrt(iterations),'--b');
yline(ax3,-2.57/sqrt(iterations),'--b');
yline(ax4,-2.57/sqrt(iterations),'--b');
yline(ax5,-2.57/sqrt(iterations),'--b');
yline(ax6,-2.57/sqrt(iterations),'--b');

mag=hypot(data_full(11,:,1),data_full(12,:,1));
dir=atan2d(data_full(12,:,1),data_full(11,:,1));

figure
plot(1:steps,current.set,1:steps,mag)
axis([0 steps 0 3])

figure
plot(1:steps,current.drift,1:steps,dir)

% std_dev1=sqrt(P_plus_out(1,:));
% std_dev2=sqrt(P_plus_out(2,:));
% std_dev3=sqrt(P_plus_out(3,:));
% std_dev4=sqrt(P_plus_out(4,:));
% std_dev5=sqrt(P_plus_out(5,:));
% std_dev6=sqrt(P_plus_out(6,:));
% 
% velocity_e=(track.stw.*cosd(track.hdg)).';
% velocity_n=(track.stw.*sind(track.hdg)).';
% current_e=(current.set.*cosd(current.drift)).';
% current_n=(current.set.*sind(current.drift)).';
% 
% figure
% plot(1:steps,data(15,:),1:steps,data(7,:),1:steps,data(15,:)+std_dev1,'--',1:steps,data(15,:)-std_dev1,'--')
% legend('Actual Easting Position','Filter Easting Position', 'Easting Position Upper Bound','Easting Position Lower Bound')
% figure
% plot(1:steps,velocity_e,1:steps,data(8,:),1:steps,velocity_e+std_dev2,'--',1:steps,velocity_e-std_dev2,'--')
% legend('Actual Easting Velocity','Filter Easting Velocity', 'Easting Velocity Upper Bound','Easting Velocity Lower Bound')
% figure
% plot(1:steps,data(16,:),1:steps,data(9,:),1:steps,data(16,:)+std_dev3,'--',1:steps,data(16,:)-std_dev3,'--')
% legend('Actual Easting Position','Filter Easting Position', 'Easting Position Upper Bound','Easting Position Lower Bound')
% figure
% plot(1:steps,velocity_n,1:steps,data(10,:),1:steps,velocity_n+std_dev4,'--',1:steps,velocity_n-std_dev4,'--')
% legend('Actual Northing Velocity','Filter Northing Velocity', 'Northing Velocity Upper Bound','Northing Velocity Lower Bound')
% figure
% plot(1:steps,current_e,1:steps,data(11,:),1:steps,current_e+std_dev5,'--',1:steps,current_e-std_dev5,'--')
% legend('Actual Easting Current','Filter Easting Current', 'Easting Current Upper Bound','Easting Current Lower Bound')
% figure
% plot(1:steps,current_n,1:steps,data(12,:),1:steps,current_n+std_dev6,'--',1:steps,current_n-std_dev6,'--')
% legend('Actual Northing Current','Filter Northing Current', 'Northing Current Upper Bound','Northing Current Lower Bound')
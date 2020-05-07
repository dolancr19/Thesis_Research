% clear variables
% clc
% 
% load 'D:\Documents\Thesis_Research\Bigelow_12_12\processed_data\1576176163.79_NAVDR.mat';
% 
% order=4;
% window=5;
% p=1;
% x=NAVDR.NAV_HEADING+360;
% dt=NAVDR.TIME(2)-NAVDR.TIME(1);
% 
% 
% dpxdtp = ddt_sgolay(x,order,window,dt,p);

% %Compute R, the covariance matrix associated with measurement error
% var_Rr=.1^2; % range
% var_Ra=.2^2; % azimuth
% var_Rs=.3^2; % speed through the water
% var_Ra1=.4^2; % heading
% R=diag([var_Rr,var_Ra,var_Rs,var_Ra1]);
% 
% v_k=randn(4,1);
% [V_R,D_R]=eig(R);
% % [~,ind] = sort(diag(D_R),'descend');
% % D_Rs = D_R(ind,ind);
% % V_Rs = V_R(:,ind);
% d=diag(D_R);
% v_k=V_R*(d.*v_k);

% set_lower=0;
% set_upper=.5;
% drift_lower=-180;
% drift_upper=180;
% 
% set=(set_upper-set_lower)*rand(1,1)+set_lower
% drift=(drift_upper-drift_lower)*rand(1,1)+drift_lower

% std_dev1=sqrt(P_plus_out(1,:));
% std_dev2=sqrt(P_plus_out(2,:));
% std_dev3=sqrt(P_plus_out(3,:));
% std_dev4=sqrt(P_plus_out(4,:));
% std_dev5=sqrt(P_plus_out(5,:));
% std_dev6=sqrt(P_plus_out(6,:));
% 
% 
% figure
% subplot(6,1,1)
% plot(1:steps,data(15,:),1:steps,data(7,:),1:steps,data(7,:)+std_dev1,'--',1:steps,data(7,:)-std_dev1,'--')
% legend('Actual Easting Position','Filter Easting Position', 'Filter Easting Position Upper Bound','Filter Easting Position Lower Bound')
% subplot(6,1,2)
% plot(1:steps,track.stw.*cosd(track.hdg),1:steps,data(8,:),1:steps,data(8,:)+std_dev2,'--',1:steps,data(8,:)-std_dev2,'--')
% legend('Actual Northing Position','Filter Northing Position', 'Filter Northing Position Upper Bound','Filter Northing Position Lower Bound')
% subplot(6,1,3)
% plot(1:steps,data(16,:),1:steps,data(9,:),1:steps,data(9,:)+std_dev3,'--',1:steps,data(9,:)-std_dev3,'--')
% legend('Actual Easting Velocity','Filter Easting Velocity', 'Filter Easting Velocity Upper Bound','Filter Easting Velocity Lower Bound')
% subplot(6,1,4)
% plot(1:steps,track.stw.*sind(track.hdg),1:steps,data(10,:),1:steps,data(10,:)+std_dev4,'--',1:steps,data(10,:)-std_dev4,'--')
% legend('Actual Northing Velocity','Filter Northing Velocity', 'Filter Northing Velocity Upper Bound','Filter Northing Velocity Lower Bound')
% subplot(6,1,5)
% plot(1:steps,current.set.*cosd(current.drift),1:steps,data(11,:),1:steps,data(11,:)+std_dev5,'--',1:steps,data(11,:)-std_dev5,'--')
% legend('Actual Easting Current','Filter Easting Current', 'Filter Easting Current Upper Bound','Filter Easting Current Lower Bound')
% subplot(6,1,6)
% plot(1:steps,current.set.*sind(current.drift),1:steps,data(12,:),1:steps,data(12,:)+std_dev6,'--',1:steps,data(12,:)-std_dev6,'--')
% legend('Actual Northing Current','Filter Northing Current', 'Filter Northing Current Upper Bound','Filter Northing Current Lower Bound')

% current_north_interp_sel=interp1(Data.BottomTrack_Time_UNIX,Data.Comp_VelNorth,local_t(sel));
% current_east_interp_sel=interp1(Data.BottomTrack_Time_UNIX,Data.Comp_VelEast,local_t(sel));
% avg_current_east_interp_sel=zeros(length(current_east_interp_sel),1);
% avg_current_north_interp_sel=zeros(length(current_east_interp_sel),1);
% for ii=1:length(current_east_interp_sel)
%     avg_current_east_interp_sel(ii)=mean(current_east_interp_sel(ii,:));
%     avg_current_north_interp_sel(ii)=mean(current_north_interp_sel(ii,:));
% end
% 
% jetyak_lat_interp_sel=interp1(Data.jetyak_Time,Data.jetyak_Lat,local_t(sel));
% jetyak_lon_interp_sel=interp1(Data.jetyak_Time,Data.jetyak_Lon,local_t(sel));
% 
% quiver(jetyak_lon_interp_sel,jetyak_lat_interp_sel, avg_current_east_interp_sel,avg_current_north_interp_sel,'AutoScaleFactor',2)

% total_time=PF_RB.Time(end)-PF_RB.Time(1);
% tenth_second=0:.1:total_time;
% tenth_tb=tenth_second+PF_RB.Time(1);

% for aa=1:length(Interp.e_gps)
%     e_source=Interp.e_gps(aa)-Interp.Range(aa)*cosd(Interp.Bearing(aa));
%     n_source=Interp.n_gps(aa)-Interp.Range(aa)*sind(Interp.Bearing(aa));
% end
% sel = zeros(length(ACOUSTIC_STDDEV.ACOUSTIC_STDDEV_1),1);
% for ii=1:length(sel)
%     if ACOUSTIC_STDDEV.ACOUSTIC_STDDEV_1(ii) <= 15 && ACOUSTIC_STDDEV.ACOUSTIC_STDDEV_2(ii) <= 15
%         sel(ii)=true;
%     else
%         sel(ii)=false;
%     end
% end
% sel=logical(sel);
% 
% for jj=1:length(sel)
%     if sel(jj)==true
%         test_out(jj)= ACOUSTIC_STDDEV.ACOUSTIC_STDDEV_1(jj);
%     else
%         test_out(jj)=0;
%     end
% end
% plot(test_out)

% test_string(1)=1;
% test_string(2)=.1;
% string(test_string(1))+'_'+string(test_string(2))

% dt=1;
% var_Q=.1^2;
% var_Qc=.01;
% Qc=diag([var_Qc var_Qc]); 
% %var_Q=.0001^2;
% G=[.5*(dt^2); dt;.5*(dt^2);dt];
% Q=G*var_Q*G';
% Q=[Q zeros(4,2);zeros(2,4) Qc]

% bias_mag=.0874;
% var_Qb = 9.5*bias_mag^2 + -0.75*bias_mag + 0.12
% %     Q(8,8)=var_Qb^2;

[~,g] = sgolay(1,7);
dt=1;
p=1;
e_vel_PF = zeros(length(e_pos_PF),1);
n_vel_PF = zeros(length(n_pos_PF),1);
e_vel_PF(:,1) = conv(e_pos_PF, factorial(p)/(-dt)^p * g(:,p+1), 'same');
n_vel_PF(:,1) = conv(n_pos_PF, factorial(p)/(-dt)^p * g(:,p+1), 'same');

e_vel_PF=e_vel_PF(4:end-3);
n_vel_PF=n_vel_PF(4:end-3);


speed_PF=hypot(e_vel_PF,n_vel_PF);
plot(speed_PF)



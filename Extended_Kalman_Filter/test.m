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

% [~,g] = sgolay(1,3);
% dt=1;
% p=1;
% e_vel_PF = zeros(length(e_pos_PF),1);
% n_vel_PF = zeros(length(n_pos_PF),1);
% e_vel_PF(:,1) = conv(e_pos_PF, factorial(p)/(-dt)^p * g(:,p+1), 'same');
% n_vel_PF(:,1) = conv(n_pos_PF, factorial(p)/(-dt)^p * g(:,p+1), 'same');
% 
% e_vel_PF=e_vel_PF(2:end-1);
% n_vel_PF=n_vel_PF(2:end-1);
% 
% speed_PF=hypot(e_vel_PF,n_vel_PF);
% 
% figure
% plot(speed_PF)
% 
% e_vel_KF = zeros(length(data(8,:)),1);
% n_vel_KF = zeros(length(data(10,:)),1);
% e_vel_KF(:,1) = conv(data(8,:), factorial(p)/(-dt)^p * g(:,p+1), 'same');
% n_vel_KF(:,1) = conv(data(10,:), factorial(p)/(-dt)^p * g(:,p+1), 'same');
% 
% e_vel_KF=e_vel_KF(2:end-1);
% n_vel_KF=n_vel_KF(2:end-1);
% 
% speed_KF=hypot(e_vel_KF,n_vel_KF);
% 
% figure
% plot(speed_KF)

% steps=3600;
% x_error=data_full(15,:,1)-data_full(7,:,1);
% y_error=data_full(16,:,1)-data_full(9,:,1); 
% figure
% plot(1:steps,x_error,1:steps,y_error)

% LBL_XY(:,1)=time_lbl;
% LBL_XY(:,2)=xs_lbl(:,1);
% LBL_XY(:,3)=ys_lbl(:,1);
% trim=zeros(length(time_lbl),1);
% for ii=1:length(trim)
%     if isnan(LBL_XY(ii,1))
%         trim(ii)=0;
%     else
%         trim(ii)=1;
%     end
% end
% trim=logical(trim);
% LBL_XY=LBL_XY(trim,:);
% LBL_XY=array2table(LBL_XY,'VariableNames',{'time','xs','ys'});
% 
% SOURCE_XY(:,1)=source_x;
% SOURCE_XY(:,2)=source_y;
% % SOURCE_XY=repmat(SOURCE_XY,length(time),1);
% SOURCE_XY=array2table(SOURCE_XY,'VariableNames',{'source_x','source_y'});
% 
% MLE_RB(:,1)=time;
% MLE_RB(:,2)=range;
% MLE_RB(:,3)=bearing;
% MLE_RB=array2table(MLE_RB,'VariableNames',{'time','range','bearing'});

% figure
% plot(piUSBL.x_absolute,piUSBL.y_absolute)
% hold on
% plot(local_x,local_y)
% % hold on
% % plot(xs_lbl(:,1),xs_lbl(:,2))

% for jj=1:length(NAVDR.NAV_HEADING)
%     if NAVDR.NAV_HEADING(jj)>180
%         NAVDR.NAV_HEADING(jj)=NAVDR.NAV_HEADING(jj)-360;
%     elseif NAVDR.NAV_HEADING(jj)<-180
%         NAVDR.NAV_HEADING(jj)=NAVDR.NAV_HEADING(jj)+360;
%     end
% end
% source_x_calc=zeros(length(Interp.Time),1);
% source_y_calc=zeros(length(Interp.Time),1);
% 
% first=1;%2559
% last=length(MLE_RB.time);%2824
% Interp.Time=MLE_RB.time(first:last);
% 
% Interp.e_gps=SOURCE_XY.source_x(first:last);
% Interp.n_gps=SOURCE_XY.source_y(first:last);
% Interp.NAV_HEADING=interp1(NAVDR.TIME,NAVDR.NAV_HEADING,Interp.Time);
% Interp.NAV_SPEED=interp1(NAVDR.TIME,NAVDR.NAV_SPEED,Interp.Time);
% Interp.Range=MLE_RB.range(first:last);
% Interp.Bearing=MLE_RB.bearing(first:last);
% 
% for ii=1:length(Interp.Time)
%     acoustic_bearing=90-(Interp.NAV_HEADING(ii)-rad2deg(Interp.Bearing(ii)));
%     if acoustic_bearing>180
%         acoustic_bearing=acoustic_bearing-360;
%     elseif acoustic_bearing<-180
%         acoustic_bearing=acoustic_bearing+360;
%     end
%     
%     source_x_calc(ii)=Interp.Range(ii)*cosd(acoustic_bearing);
%     source_y_calc(ii)=Interp.Range(ii)*sind(acoustic_bearing);
% end

% trim_z=zeros(1,length(z_k_out));
% for ii=1:length(z_k_out)
%     if z_k_out(1,ii)==0
%         trim_z(ii)=0;
%     else
%         trim_z(ii)=1;
%     end
% end
% trim_z=logical(trim_z);
% z_k_trim=z_k_out(:,trim_z);

% R=[.1 .2;.3 .4];
% R_inv=inv(R)

% x_k_minus=[2;0;2;0;0;0;1];
% %Calculate measurement vector, both linear and non-linear elements
%         h=[x_k_minus(7,1)*hypot(x_k_minus(2,1),x_k_minus(4,1));atan2d(x_k_minus(4,1),x_k_minus(2,1))];
% 
%         %Calculate measurement mapping matrix
%         if hypot(x_k_minus(2,1),x_k_minus(4,1))==0
%             stw_denom=0;
%             hdg_dx=0;
%             hdg_dy=0;
%         else
%             stw_denom=x_k_minus(7,1)/hypot(x_k_minus(2,1),x_k_minus(4,1));
%             stw_dx=x_k_minus(2,1)*stw_denom;
%             stw_dy=x_k_minus(4,1)*stw_denom;
% 
%             hdg_dx=-1*x_k_minus(4,1)/(x_k_minus(2,1)^2+x_k_minus(4,1)^2);
%             hdg_dy=1/(x_k_minus(2,1)+(x_k_minus(4,1)^2/x_k_minus(2,1)));
%         end
% 
%         H=[0 stw_dx 0 stw_dy 0 0 1;0 hdg_dx 0 hdg_dy 0 0 0];
% acoustic_count=1;
% ii=5553;
% while NAVDR.TIME(ii)-LBL_NAV.time(acoustic_count)>0
%     acoustic_count=acoustic_count+1;
% end
% acoustic_count=acoustic_count-1
% NAVDR.TIME(ii)-MLE_RB.time(acoustic_count)

% LBL_NAV(:,1)=time_lbl;
% LBL_NAV(:,2)=xs_lbl(:,1);
% LBL_NAV(:,3)=ys_lbl(:,1);
% LBL_NAV=array2table(LBL_NAV,'VariableNames',{'time','xs','ys'});

% MLE_BRG_deg=zeros(length(MLE_RB.time),1);
% for ii=1:length(MLE_RB.time)
%     acoustic_bearing=rad2deg(MLE_RB.bearing(ii));
%     if acoustic_bearing>180
%         acoustic_bearing=acoustic_bearing-360;
%     elseif acoustic_bearing<-180
%         acoustic_bearing=acoustic_bearing+360;
%     end
%     
%     MLE_BRG_deg(ii)=acoustic_bearing;
%     
% end
% ll=1;
% for mm=1:length(data(8,:))
%     if mod(mm,50)==0
%         quiver_interp(1,ll)=data(8,mm);
%         quiver_interp(2,ll)=data(10,mm);
%         quiver_interp(3,ll)=data(12,mm);
%         quiver_interp(4,ll)=data(13,mm);
%         ll=ll+1;
%     end
% end
% 
% scale_factor=30;
% figure(11)
% reset(gca)
% 
% quiver(quiver_interp(1,:),quiver_interp(2,:),scale_factor*quiver_interp(3,:),scale_factor*quiver_interp(4,:),'AutoScale','off')
% % hold on
% % quiver(Interp.e_gps,Interp.n_gps,scale_factor*Interp.avg_VelEast,scale_factor*Interp.avg_VelNorth,'AutoScale','off')
% hold on
% quiver(-50,-160,scale_factor*0,scale_factor*1,'AutoScale','off')
% hold on
% plot(data(8,first:end),data(10,first:end))
% axis equal
% % legend('Calculated current','Measured current','1 m/s')
% legend('Calculated current','1 m/s', 'Filter calculated track')
% xlabel('Easting position (m)')
% ylabel('Northing position (m)')

directory="D:\Documents\Thesis_Research\MIT_Sailing\processed_data\quokka\";
start_time=1536932173.36;
%% Define origin of local coordinate system
lat=42.35846;
lon=-71.08759;

%% Define UTM structure
zone=utmzone(lat,lon);
utmstruct = defaultm('utm'); 
utmstruct.zone = zone;  
utmstruct.geoid = wgs84Ellipsoid;
utmstruct = defaultm(utmstruct);

%% Convert origin to UTM
[x,y] = mfwdtran(utmstruct,lat,lon);
% MLE_utm=[z_k_trim(1,:)+x;z_k_trim(2,:)+y];
% LBL_utm=[LBL_NAV.xs(569:3810)+x LBL_NAV.ys(569:3810)+y]; 
% LBL_utm=[LBL_XY.xs(569:3810)+x LBL_XY.ys(569:3810)+y]; 
source_utm=[SOURCE_XY.source_x+x SOURCE_XY.source_y+y];


%% Convert UTM to LL
% [lat1,lon1] = minvtran(utmstruct,MLE_utm(1,:),MLE_utm(2,:));
% [lat1,lon1] = minvtran(utmstruct,LBL_utm(:,1),LBL_utm(:,2));
[lat1,lon1] = minvtran(utmstruct,source_utm(:,1),source_utm(:,2));
s=geoshape(lat1,lon1);
filename_kml=directory + string(start_time) + "_source.kml";
% filename_kml=directory + string(start_time) + "_LBL1.kml";
kmlwrite(filename_kml,s);
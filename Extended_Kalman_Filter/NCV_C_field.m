%function [data,epsilon_v]= NCV_field(NAVDR, SOURCEXY, ACOUSTICRB, e_gps_interp, n_gps_interp) %MOOS data
function [data,epsilon_v,P_minus_out,H_out,S_out,K_out,z_k_out,h_out,F_out]=NCV_C_field(Interp) %PF data
steps=length(Interp.Time);
data=zeros(23,steps); %matrix for post-processing data

%Define boundary conditions
acoustic_bearing=90-(Interp.NAV_HEADING(1)-Interp.Bearing(1));
if acoustic_bearing>180
    acoustic_bearing=acoustic_bearing-360;
elseif acoustic_bearing<-180
    acoustic_bearing=acoustic_bearing+360;
end
source_x=Interp.Range(1)*cosd(acoustic_bearing);
source_y=Interp.Range(1)*sind(acoustic_bearing);

e0_pos=Interp.e_gps(1)-source_x;
n0_pos=Interp.n_gps(1)-source_y;

stw0=Interp.NAV_SPEED(1);
hdg0=90-Interp.NAV_HEADING(1);

e0_vel=stw0*cosd(hdg0);
n0_vel=stw0*sind(hdg0);

%r_k=zeros(4,1); % initialize innovation vector
%S=zeros(length(r_k)); % initialize innovation covariance matrix

%Initialize state vectors
x_k_minus=[e0_pos;e0_vel;n0_pos;n0_vel;0;0];
data(1:6,1)=x_k_minus;
x_k_plus=x_k_minus;
data(7:12,1)=x_k_plus;
z_k_out=zeros(4,steps);
z_k_out(:,1)=x_k_minus(1:4,1);
%Define Q, the covariance matrix for noise associated with the state vector
var_Qp=5; % position
var_Qw=.1; % water referenced velocity
var_Qc=.01; % current velocity
Q=diag([var_Qp^2,var_Qw^2,var_Qp^2,var_Qw^2,var_Qc^2,var_Qc^2]);

%Compute R, the covariance matrix associated with measurement error
var_Rr_r=10^2; % range
var_Ra_r=deg2rad(1^2); % azimuth
var_Rs_r=.1^2; % speed through the water
var_Ra1_r=deg2rad(1^2); % heading
%R=diag([var_Rr,var_Ra,var_Rs,var_Ra1]);

%Initialize error covariance matrices
var_Pr=10^2; % position
var_Ps=10^2; % velocity
P_minus=diag([var_Pr,var_Ps,var_Pr,var_Ps,var_Ps,var_Ps]);
P_plus=P_minus;

H=[1 0 0 0 0 0;0 1 0 0 0 0;0 0 1 0 0 0;0 0 0 1 0 0];

%Initialize vectors for consistency checks
%epsilon=zeros(1,steps);
epsilon_v=zeros(1,steps);
%mu=zeros(length(x_k_plus),steps);
K_out=zeros(6,4,steps);
H_out=zeros(4,6,steps);
S_out=zeros(4,4,steps);
P_minus_out=zeros(6,6,steps);

h_out=zeros(4,steps);
F_out=zeros(6,6,steps);
%% Execute filter
for ii=2:steps
    dt=Interp.Time(ii)-Interp.Time(ii-1);
    F=[1 dt 0 0 dt 0;0 1 0 0 0 0;0 0 1 dt 0 dt;0 0 0 1 0 0;0 0 0 0 1 0;0 0 0 0 0 1]; % state transition matrix
    F_out(:,:,ii)=F;
    %Calculate state estimate and error covariance matrix for next step
    x_k_minus=F*x_k_plus;
    P_minus=F*P_plus*F.'+Q;
    P_minus_out(:,:,ii)=P_minus;    
    data(1:6,ii)=x_k_minus;
    
    acoustic_bearing=90-(Interp.NAV_HEADING(ii)-Interp.Bearing(ii));
    if acoustic_bearing>180
        acoustic_bearing=acoustic_bearing-360;
    elseif acoustic_bearing<-180
        acoustic_bearing=acoustic_bearing+360;
    end
    
    range=Interp.Range(ii);
    azi_r=deg2rad(acoustic_bearing);
    stw=Interp.NAV_SPEED(ii);
    hdg_r=deg2rad(90-Interp.NAV_HEADING(ii));
    
    
    mu_t=[range*cos(azi_r)*(exp(-1*var_Ra_r)-exp(-.5*var_Ra_r));stw*cos(hdg_r)*(exp(-1*var_Ra1_r)-exp(-.5*var_Ra1_r));range*sin(azi_r)*(exp(-1*var_Ra_r)-exp(-.5*var_Ra_r));stw*sin(hdg_r)*(exp(-1*var_Ra1_r)-exp(-.5*var_Ra1_r))];
    
    R(1,1)=(range^2*exp(-2*var_Ra_r)*((cos(azi_r)^2)*(cosh(2*var_Ra_r)-cosh(var_Ra_r))+(sin(azi_r)^2*(sinh(2*var_Ra_r)-sinh(var_Ra_r)))))+(var_Rr_r*exp(-2*var_Ra_r)*((cos(azi_r)^2)*(2*cosh(2*var_Ra_r)-cosh(var_Ra_r))+(sin(azi_r)^2*(2*sinh(2*var_Ra_r)-sinh(var_Ra_r)))));
    R(2,2)=(stw^2*exp(-2*var_Ra1_r)*((cos(hdg_r)^2)*(cosh(2*var_Ra1_r)-cosh(var_Ra1_r))+(sin(hdg_r)^2*(sinh(2*var_Ra1_r)-sinh(var_Ra1_r)))))+(var_Rs_r*exp(-2*var_Ra1_r)*((cos(hdg_r)^2)*(2*cosh(2*var_Ra1_r)-cosh(var_Ra1_r))+(sin(hdg_r)^2*(2*sinh(2*var_Ra1_r)-sinh(var_Ra1_r)))));
    R(3,3)=(range^2*exp(-2*var_Ra_r)*((sin(azi_r)^2)*(cosh(2*var_Ra_r)-cosh(var_Ra_r))+(cos(azi_r)^2*(sinh(2*var_Ra_r)-sinh(var_Ra_r)))))+(var_Rr_r*exp(-2*var_Ra_r)*((sin(azi_r)^2)*(2*cosh(2*var_Ra_r)-cosh(var_Ra_r))+(cos(azi_r)^2*(2*sinh(2*var_Ra_r)-sinh(var_Ra_r)))));
    R(4,4)=(stw^2*exp(-2*var_Ra1_r)*((sin(hdg_r)^2)*(cosh(2*var_Ra1_r)-cosh(var_Ra1_r))+(cos(hdg_r)^2*(sinh(2*var_Ra1_r)-sinh(var_Ra1_r)))))+(var_Rs_r*exp(-2*var_Ra1_r)*((sin(hdg_r)^2)*(2*cosh(2*var_Ra1_r)-cosh(var_Ra1_r))+(cos(hdg_r)^2*(2*sinh(2*var_Ra1_r)-sinh(var_Ra1_r)))));

    R(1,3)=sin(azi_r)*cos(azi_r)*exp(-4*var_Ra_r)*(var_Rr_r+(range^2+var_Rr_r)*(1-exp(var_Ra_r)));
    R(3,1)=R(1,3);
    
    R(2,4)=sin(hdg_r)*cos(hdg_r)*exp(-4*var_Ra1_r)*(var_Rs_r+(stw^2+var_Ra1_r)*(1-exp(var_Ra1_r)));
    R(4,2)=R(2,4);
    
    z_k=[Interp.e_gps(ii)-Interp.Range(ii)*cosd(acoustic_bearing);stw*cos(hdg_r);Interp.n_gps(ii)-Interp.Range(ii)*sind(acoustic_bearing);stw*sin(hdg_r)]-mu_t;
    
    %z_k=[Interp.Range(ii);acoustic_bearing;Interp.NAV_SPEED(ii);90-Interp.NAV_HEADING(ii)];
    z_k_out(:,ii)=z_k;
%     %Calculate non-linear measurement
%     h=[hypot(Interp.e_gps(ii)-x_k_minus(1,1),Interp.n_gps(ii)-x_k_minus(3,1));atan2d(Interp.n_gps(ii)-x_k_minus(3,1),Interp.e_gps(ii)-x_k_minus(1,1));hypot(x_k_minus(2,1),x_k_minus(4,1));atan2d(x_k_minus(4,1),x_k_minus(2,1))];
%     h_out(:,ii)=h;
%     %Calculate measurement mapping Jacobian 
%     range_denom=1/hypot(Interp.e_gps(ii)-x_k_minus(1,1),Interp.n_gps(ii)-x_k_minus(3,1));
%     range_dx=(Interp.e_gps(ii)-x_k_minus(1,1))*range_denom;
%     range_dy=(Interp.n_gps(ii)-x_k_minus(3,1))*range_denom;
% 
%     azi_dx=-1*(Interp.n_gps(ii)-x_k_minus(3,1))/((Interp.e_gps(ii)-x_k_minus(1,1))^2+(Interp.n_gps(ii)-x_k_minus(3,1))^2);
%     azi_dy=1/((Interp.e_gps(ii)-x_k_minus(1,1))+((Interp.n_gps(ii)-x_k_minus(3,1))^2/(Interp.e_gps(ii)-x_k_minus(1,1))));
% 
%     stw_denom=1/hypot(x_k_minus(2,1),x_k_minus(4,1));
%     stw_dx=x_k_minus(2,1)*stw_denom;
%     stw_dy=x_k_minus(4,1)*stw_denom;
% 
%     hdg_dx=-1*x_k_minus(4,1)/(x_k_minus(2,1)^2+x_k_minus(4,1)^2);
%     hdg_dy=1/(x_k_minus(2,1)+(x_k_minus(4,1)^2/x_k_minus(2,1)));
% 
%     H=[range_dx 0 range_dy 0 0 0;azi_dx 0 azi_dy 0 0 0;0 stw_dx 0 stw_dy 0 0;0 hdg_dx 0 hdg_dy 0 0];
%     H_out(:,:,ii)=H;
    %Calculate innovation
    r_k=z_k-H*x_k_minus;
%     r_k(1,1)=z_k(1,1)-h(1,1);
%     r_k(2,1)=angdiff(deg2rad(h(2,1)),deg2rad(z_k(2,1)));
%     r_k(3,1)=z_k(3,1)-h(3,1);
%     r_k(4,1)=angdiff(deg2rad(h(4,1)),deg2rad(z_k(4,1)));

    %Calculate innovation covariance
    S=H*P_minus*H.'+R;
    S_out(:,:,ii)=S;
    %Update Kalman Gain after estimate vector calculation
    K=(P_minus*H.')/S;
    K_out(:,:,ii)=K;

    %Update state vector
    x_k_plus=x_k_minus+K*r_k;

    %Update error covariance matrix with new Kalman Gain
    %P_plus=(eye(4)-K*H)*P_minus;

    %Joseph stabilized equation
    P_plus=(eye(6)-K*H)*P_minus*(eye(6)-K*H)'+K*R*K';
%     else
%         x_k_plus=x_k_minus;
%         P_plus=P_minus;
%     end
%     %Calculate the normalized estimation error squared (NEES)
%     x_err=x_act-x_k_plus;
%     epsilon(ii)=x_err.'*(P_plus\x_err);
%     
%     %Calculate the normalized mean estimation error (NMEE)
%     for jj=1:length(x_k_plus)
%         mu(jj,ii)=x_err(jj,1)/sqrt(P_plus(jj,jj));
%     end
    
    %Calculate the normalized innovation squared (NIS)
    epsilon_v(ii)=r_k.'*(S\r_k);
        
    data(7:12,ii)=x_k_plus;
    data(13,ii)=hypot(x_k_plus(2,1),x_k_plus(4,1));
    data(14,ii)=atan2d(x_k_plus(4,1),x_k_plus(2,1));
    data(15,ii)=norm(P_plus);
    %data(16:19,ii)=S;
    data(20:23,ii)=r_k;
end
epsilon_v(isnan(epsilon_v))=0;
epsilon_v_bar=(1/steps)*sum(epsilon_v);

% sum_kj=zeros(length(r_k),1);
% sum_kk=zeros(length(r_k),1);
% sum_jj=zeros(length(r_k),1);
% for jj=1:steps-1
%     sum_kj(1)=sum_kj(1)+(data(20,jj)*data(20,jj+1));
%     sum_kj(2)=sum_kj(2)+(data(21,jj)*data(21,jj+1));
%     sum_kj(3)=sum_kj(3)+(data(22,jj)*data(22,jj+1));
%     sum_kj(4)=sum_kj(4)+(data(23,jj)*data(23,jj+1));
%     sum_kk(1)=sum_kk(1)+(data(20,jj))^2;
%     sum_kk(2)=sum_kk(2)+(data(21,jj))^2;
%     sum_kk(3)=sum_kk(3)+(data(22,jj))^2;
%     sum_kk(4)=sum_kk(4)+(data(23,jj))^2;
%     sum_jj(1)=sum_jj(1)+(data(20,jj+1))^2;
%     sum_jj(2)=sum_jj(2)+(data(20:23,jj+1))^2;
%     sum_jj=sum_jj+(data(20:23,jj+1))^2;
%     sum_jj=sum_jj+(data(20:23,jj+1))^2;
% end
% rho_bar=zeros(4,1);
% rho_bar(1,1)=sum_kj(1)*sqrt(sum_kk(1)*sum_jj(1));
% rho_bar(2,1)=sum_kj(2)*sqrt(sum_kk(2)*sum_jj(2));
% rho_bar(3,1)=sum_kj(3)*sqrt(sum_kk(3)*sum_jj(3));
% rho_bar(4,1)=sum_kj(4)*sqrt(sum_kk(4)*sum_jj(4));

%% Plot results

mag=hypot(data(11,:),data(12,:));
dir=atan2d(data(12,:),data(11,:));

% e_source=zeros(1,steps);
% n_source=zeros(1,steps);
% 
% for aa=1:length(Interp.e_gps)
%     e_source(aa)=Interp.e_gps(aa)-Interp.Range(aa)*cosd(Interp.Bearing(aa));
%     n_source(aa)=Interp.n_gps(aa)-Interp.Range(aa)*sind(Interp.Bearing(aa));
% end

figure(1)
plot(1:steps,data(15,1:steps))
title('Error Covariance Matrix Norm vs Filter Step')
xlabel('Filter Step')
ylabel('Error Covariance Matrix Norm')

% figure
% plot(1:steps,data(18,:))
% %axis([0 steps -5 5])
% title('Filter Error vs Filter Step')
% xlabel('Filter Step')
% ylabel('Filter Error (m)')
% 
% figure
% plot(1:steps,data(19,:))
% %axis([0 steps -5 5])
% title('Filter Error vs Filter Step')
% xlabel('Filter Step')
% ylabel('Filter Error (degrees)')
% 
%%
figure(2)
plot(data(7,:),data(9,:),z_k_out(1,:),z_k_out(3,:),Interp.e_gps,Interp.n_gps)
axis equal
title('Filter Output Position')
xlabel('x position (m)')
ylabel('y position (m)')
legend('Filter Output','Input Measurement','Jetyak Position')
%%
figure(3)
reset(gca)
yyaxis left
plot(1:steps,data(20,:),1:steps,data(22,:))
title('Innovations')
xlabel('Step')
ylabel('Position (m)')
yyaxis right
plot(1:steps,data(21,:),1:steps,data(23,:))
ylabel('Angle (degrees)')
legend('Range','Azimuth','STW','Heading')

figure(4)
plot(1:steps,mag)
% axis([0 steps 0 3])
% 
figure(5)
plot(1:steps,dir)

% figure
% plot(epsilon)
% 
% figure
% ax1 = subplot(6,1,1);
% x = 1:steps;
% y1 = mu(1,:); 
% plot(ax1,x,y1)
% axis([0 steps -2.25 2.25])
% 
% ax2 = subplot(6,1,2);
% x = 1:steps;
% y2 = mu(2,:); 
% plot(ax2,x,y2)
% axis([0 steps -2.25 2.25])
% 
% ax3 = subplot(6,1,3);
% x = 1:steps;
% y3 = mu(3,:); 
% plot(ax3,x,y3)
% axis([0 steps -2.25 2.25])
% 
% ax4 = subplot(6,1,4);
% x = 1:steps;
% y4 = mu(4,:); 
% plot(ax4,x,y4)
% axis([0 steps -2.25 2.25])
% 
% ax5 = subplot(6,1,5);
% x = 1:steps;
% y5 = mu(5,:); 
% plot(ax5,x,y5)
% axis([0 steps -2.25 2.25])
% 
% ax6 = subplot(6,1,6);
% x = 1:steps;
% y6 = mu(6,:); 
% plot(ax6,x,y6)
% axis([0 steps -2.25 2.25])
% 
% yline(ax1,1.96,'--b');
% yline(ax2,1.96,'--b');
% yline(ax3,1.96,'--b');
% yline(ax4,1.96,'--b');
% yline(ax5,1.96,'--b');
% yline(ax6,1.96,'--b');
% yline(ax1,-1.96,'--b');
% yline(ax2,-1.96,'--b');
% yline(ax3,-1.96,'--b');
% yline(ax4,-1.96,'--b');
% yline(ax5,-1.96,'--b');
% yline(ax6,-1.96,'--b');
% 
figure(6)
plot(epsilon_v_bar)


 %% Plot quiver
reset(gca);
figure(7)
quiver(data(7,:),data(9,:),data(11,:),data(12,:))
% figure
% plot(1:steps,data(14,:),1:steps,track(:,2))

end
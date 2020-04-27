function [data,epsilon,epsilon_v,mu, K_out, H_out, P_plus_out, P_minus_out]= NCV_bias(freq, steps, x0_pos, y0_pos, stw0, hdg0, track, current,track0)
dt=1/freq; %Define time step size
data=zeros(29,steps); %matrix for post-processing data

x0_vel=stw0*cosd(hdg0);
y0_vel=stw0*sind(hdg0);

r_k=zeros(4,1); % initialize innovation vector
%S=zeros(length(r_k)); % initialize innovation covariance matrix

%Initialize state vectors
x_k_minus=[x0_pos;x0_vel;y0_pos;y0_vel;0;0;1]; 
x_k_plus=x_k_minus;

%Define Q, the covariance matrix for noise associated with the state vector
var_Qp=.575; % position
var_Qw=.035; % water referenced velocity
var_Qc=.01; % current velocity
var_Qb=.001; % bias coefficient
Q=diag([var_Qp^2,var_Qw^2,var_Qp^2,var_Qw^2,var_Qc^2,var_Qc^2,var_Qb^2]);

%Compute R, the covariance matrix associated with measurement error
var_Rr=10^2; % range
var_Ra=1^2; % azimuth
var_Rs=.1^2; % speed through the water
var_Ra1=1^2; % heading
R=diag([var_Rr,var_Ra,var_Rs,var_Ra1]);

%Initialize error covariance matrices
var_Pr=1^2; % position
var_Ps=.1^2; % velocity
var_Pb=1^2; % bias
P_minus=diag([var_Pr,var_Ps,var_Pr,var_Ps,var_Ps,var_Ps,var_Pb]);
P_plus=P_minus;

%Initialize measurement vector
range_x=0;
range_y=0;
% range=0;
% azi=0;
% stw=0;
% hdg=0;
%z_k=[range;azi;stw;hdg];



F=[1 dt 0 0 dt 0 0;0 1 0 0 0 0 0;0 0 1 dt 0 dt 0;0 0 0 1 0 0 0;0 0 0 0 1 0 0;0 0 0 0 0 1 0;0 0 0 0 0 0 1]; % state transition matrix

%Initialize vectors for consistency checks
epsilon=zeros(1,steps);
epsilon_v=zeros(1,steps);
mu=zeros(length(x_k_plus),steps);
K_out=zeros(length(x_k_plus),steps);
H_out=zeros(4,length(x_k_plus),steps);
P_minus_out=zeros(length(x_k_plus),steps);
P_plus_out=zeros(length(x_k_plus),steps);
%% Execute filter
for ii=1:steps
    %Calculate state estimate and error covariance matrix for next step
    x_k_minus=F*x_k_plus;
    bias_noise=randn(1,1)*var_Qb;
    x_k_minus(7,1)=x_k_minus(7,1)+bias_noise;
    P_minus=F*P_plus*F.'+Q;
    P_minus_out(:,ii)=diag(P_minus);    
    data(1:7,ii)=x_k_minus;
               
    %Calculate noisy measurement vector
    x_vel=track0.stw(ii)*cosd(track0.hdg(ii))+current.set(ii)*cosd(current.drift(ii));
    y_vel=track0.stw(ii)*sind(track0.hdg(ii))+current.set(ii)*sind(current.drift(ii));
    range_x=range_x+x_vel*dt;
    range_y=range_y+y_vel*dt;
    range=hypot(range_x,range_y);
    azi=atan2d(range_y,range_x);
    
    x_act=[range_x;track.stw(ii)*cosd(track.hdg(ii));range_y;track.stw(ii)*sind(track.hdg(ii));current.set(ii)*cosd(current.drift(ii));current.set(ii)*sind(current.drift(ii));1];
    
    %Executes the filter at the specified frequency with one second
    %measurement updates.  Simulates a 1 Hz iUSBL transmission.
%    if mod(ii,freq)==0
        %Calculate noisy measurement vector
        v_k=randn(4,1);
        [V_R,D_R]=eig(R);
%         [~,ind] = sort(diag(D_R),'descend');
%         D_Rs = D_R(ind,ind);
%         V_Rs = V_R(:,ind);
        d=diag(D_R);
        d=sqrt(d);
        v_k=V_R*(d.*v_k);
        z_k=[range;azi;track.stw(ii);track.hdg(ii)]+v_k;

        %Calculate non-linear measurement
        h=[hypot(x_k_minus(1,1),x_k_minus(3,1));atan2d(x_k_minus(3,1),x_k_minus(1,1));x_k_minus(7,1)*hypot(x_k_minus(2,1),x_k_minus(4,1));atan2d(x_k_minus(4,1),x_k_minus(2,1))];

        %Calculate measurement mapping Jacobian 
        range_denom=1/hypot(x_k_minus(1,1),x_k_minus(3,1));
        range_dx=x_k_minus(1,1)*range_denom;
        range_dy=x_k_minus(3,1)*range_denom;

        azi_dx=-1*x_k_minus(3,1)/(x_k_minus(1,1)^2+x_k_minus(3,1)^2);
        azi_dy=1/(x_k_minus(1,1)+(x_k_minus(3,1)^2/x_k_minus(1,1)));

        stw_denom=x_k_minus(7,1)/hypot(x_k_minus(2,1),x_k_minus(4,1));
        stw_dx=x_k_minus(2,1)*stw_denom;
        stw_dy=x_k_minus(4,1)*stw_denom;

        hdg_dx=-1*x_k_minus(4,1)/(x_k_minus(2,1)^2+x_k_minus(4,1)^2);
        hdg_dy=1/(x_k_minus(2,1)+(x_k_minus(4,1)^2/x_k_minus(2,1)));

        H=[range_dx 0 range_dy 0 0 0 0;azi_dx 0 azi_dy 0 0 0 0;0 stw_dx 0 stw_dy 0 0 0;0 hdg_dx 0 hdg_dy 0 0 0];
        H_out(:,:,ii)=H;
        %Calculate innovation
%        r_k=z_k-h;
        r_k(1,1)=z_k(1,1)-h(1,1);
        r_k(2,1)=angdiff(deg2rad(h(2,1)),deg2rad(z_k(2,1)));
        r_k(3,1)=z_k(3,1)-h(3,1);
        r_k(4,1)=angdiff(deg2rad(h(4,1)),deg2rad(z_k(4,1)));

        %Calculate innovation covariance
        S=H*P_minus*H.'+R;

        %Update Kalman Gain after estimate vector calculation
        K=(P_minus*H.')/S;
        K_out(:,ii)=K*r_k;

        %Update state vector
        x_k_plus=x_k_minus+K*r_k;

        %Update error covariance matrix with new Kalman Gain
%        P_plus=(eye(6)-K*H)*P_minus;

        %Joseph stabilized equation
        P_plus=(eye(length(x_k_plus))-K*H)*P_minus*((eye(length(x_k_plus))-K*H).')+K*R*K.';
%    else
%        x_k_plus=x_k_minus;
%        P_plus=P_minus;
%    end
    P_plus_out(:,ii)=diag(P_plus);
    %Calculate the normalized estimation error squared (NEES)
    x_err=x_act-x_k_plus;
    epsilon(ii)=x_err.'*inv(P_plus)*x_err;
    
    %Calculate the normalized mean estimation error (NMEE)
    for jj=1:length(x_k_plus)
        mu(jj,ii)=x_err(jj,1)/sqrt(P_plus(jj,jj));
    end
    
    %Calculate the normalized innovation squared (NIS)
    epsilon_v(ii)=r_k.'*(S\r_k);
        
    data(8:14,ii)=x_k_plus;
    data(15,ii)=hypot(x_k_plus(2,1),x_k_plus(4,1));
    data(16,ii)=atan2d(x_k_plus(4,1),x_k_plus(2,1));
    data(17,ii)=range_x;
    data(18,ii)=range_y;
    data(19,ii)=norm(P_plus);
    data(20,ii)=range-hypot(x_k_plus(1,1),x_k_plus(3,1));
    data(21,ii)=azi-atan2d(x_k_plus(3,1),x_k_plus(1,1));
    data(22,ii)=S(1,1);
    data(23,ii)=S(2,1);
    data(24,ii)=S(3,1);
    data(25,ii)=S(4,1);
    data(26,ii)=r_k(1,1);
    data(27,ii)=r_k(2,1);
    data(28,ii)=r_k(3,1);
    data(29,ii)=r_k(4,1);
end
epsilon_v(isnan(epsilon_v))=0;
%% Plot results

% mag=hypot(data(11,:),data(12,:));
% dir=atan2d(data(12,:),data(11,:));
% 
% figure
% plot(1:steps,data(17,1:steps))
% title('Error Covariance Matrix Norm vs Filter Step')
% xlabel('Filter Step')
% ylabel('Error Covariance Matrix Norm')
% 
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
% figure
% plot(data(7,:),data(9,:),data(15,:),data(16,:))
% axis equal
% title('Acutal Position vs Filter Output Position')
% xlabel('x position (m)')
% ylabel('y position (m)')
% legend('Filter Output','Actual Position','location','southeast')
% 
% % figure
% % yyaxis left
% % plot(1:steps,data(24,:),1:steps,data(26,:))
% % yyaxis right
% % plot(1:steps,data(25,:),1:steps,data(27,:))
% 
% figure
% plot(1:steps,set,1:steps,mag)
% axis([0 steps 0 3])
% 
% figure
% plot(1:steps,drift,1:steps,dir)
% 
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
% figure
% plot(epsilon_v)
% 
% figure
% plot(1:steps,data(14,:),1:steps,track(:,2))

end
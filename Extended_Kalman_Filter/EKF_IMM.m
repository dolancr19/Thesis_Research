clear variables;
clc

%Define time step size
freq=10;
dt=1/freq;

%Define # of filter steps desired
steps=4000*freq;

track=build_track(freq,steps);
%track(:,2)=track(:,2)-90;
%Initialize data matrix
data=zeros(28,steps);
avg=zeros(4,1);

%Define boundary conditions
x0_pos=0;
y0_pos=0;
stw=.50;
hdg=30;

x0_vel=stw*cosd(hdg);
y0_vel=stw*sind(hdg);
S=[0;0;0;0];
r_k=[0;0;0;0];

%Initialize estimated state vector
x_k_minus=[x0_pos;y0_pos;x0_vel;y0_vel;0;0];
%Compute Q, the covariance matrix for noise associated with the state vector
var_Qp=.05;
var_Qc=.05;
var_Q=diag([var_Qp^2,var_Qp^2,var_Qp^2,var_Qp^2,var_Qc^2,var_Qc^2]);
% var_Q=0;
G=[.5*dt^2 0 0;.5*dt^2 0 0;0 dt 0;0 dt 0;0 0 dt;0 0 dt];
Q=G*G'*var_Q;

%Compute C, the linear transformation matrix for white noise
[U_P,T_P,V_P]=svd(Q);
S_P=sqrt(T_P);
C=U_P*S_P;

%Compute R, the covariance matrix associated with measurement error
var_Rr=8^2;
var_Ra=1^2;
var_Rs=.1^2;
var_Ra1=.5^2;
%var_Rs=.1^2;
%var_Rh=.5^2;
% var_Rr=.00001;
% var_Ra=.00001;
% var_Rs=.00001;
%var_Rh=.00001;
R=diag([var_Rr,var_Ra,var_Rs,var_Ra1]);

%Initialize error covariance matrix
var_Pr=10^2;
var_Ps=1^2;
%var_Pa=1^2;
P_minus=diag([var_Pr,var_Pr,var_Ps,var_Ps,var_Ps,var_Ps]);
P_plus=P_minus;

%Initialize measurement vector
range_x=0;
range_y=0;
range=0;
azi=0;
z_k=[range;azi];

%Initialize state vector
x_k_plus=x_k_minus;

%F=[1 0 0 0 dt 0 0 0;0 1 0 0 0 dt 0 0;0 0 0 0 1 0 0 0;0 0 0 0 0 1 0 0;0 0 0 0 1 0 0 0;0 0 0 0 0 1 0 0;0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0];
F=[1 0 dt 0 dt 0;0 1 0 dt 0 dt;0 0 1 0 0 0;0 0 0 1 0 0;0 0 0 0 1 0;0 0 0 0 0 1];

%Calculate current values
set=zeros(1,steps);
drift=zeros(1,steps);
% for ii=1:steps
%     
% end

for ii=1:steps
    
    stw=track(ii,1);
    hdg=track(ii,2);
    
    %Specify currents for simulation
%     set(ii)=.1*sin(ii*1.454e-4)+.5;
    set(ii)=0.5;
    drift(ii)=60;
    
    x_vel=stw*cosd(hdg)+set(ii)*cosd(drift(ii));
    y_vel=stw*sind(hdg)+set(ii)*sind(drift(ii));
    
    %Calculate noisy estimate vector for next step
    %u_k=[stw*dt*cosd(hdg);stw*dt*sind(hdg);stw*cosd(hdg);stw*sind(hdg);0;0];
    w_k=randn(6,1);
    w_k=C*w_k;
%     w_k(5,1)=.02;
%     w_k(6,1)=.02;
    
    %F=[1 0 0 0 dt 0;0 1 0 0 0 dt;0 0 0 0 1 0;0 0 0 0 0 1;0 0 0 0 1 0;0 0 0 0 0 1];
    %f=[x_k_plus(1,1)+(x_k_plus(5,1)*dt); x_k_plus(2,1)+(x_k_plus(6,1)*dt); x_k_plus(5,1); x_k_plus(6,1); x_k_plus(5,1); x_k_plus(6,1)];
    %F=[1 0 dt 0 0 0;0 1 0 dt 0 0;0 0 1 0 0 0;0 0 0 1 0 0;0 0 0 0 1 0;0 0 0 0 0 1];
    %x_k_minus=f+u_k+w_k;
    %x_k_minus=F*x_k_plus;%+w_k;
    %x_k_minus=F*x_k_plus+u_k+w_k;
    x_k_minus=F*x_k_plus+w_k;
    
   
    
    
    %Calculate error covariance matrix for next step
    P_minus=F*P_plus*F'+Q;
        
    data(1,ii)=x_k_minus(1,1);
    data(2,ii)=x_k_minus(2,1);
    data(3,ii)=x_k_minus(3,1);
    data(4,ii)=x_k_minus(4,1);
    data(5,ii)=x_k_minus(5,1);
    data(6,ii)=x_k_minus(6,1);
    
    
        
    %Calculate noisy measurement vector
    
    range_x=range_x+x_vel*dt;
    range_y=range_y+y_vel*dt;
    range=hypot(range_x,range_y);
    azi=atan2d(range_y,range_x);
    
    x_act=[range_x;range_y;stw*cosd(hdg);stw*sind(hdg);set(ii)*cosd(drift(ii));set(ii)*sind(drift(ii))];
    
    if mod(ii,freq)==0
    v_k=randn(4,1);
    [V_R,D_R]=eig(R);
    [~,ind] = sort(diag(D_R),'descend');
    D_Rs = D_R(ind,ind);
    V_Rs = V_R(:,ind);
    v_k=V_Rs*D_Rs*v_k;
    z_k=[range;azi;stw;hdg]+v_k;
    
    %Calculate non-linear measurement
    h=[hypot(x_k_minus(1,1),x_k_minus(2,1));atan2d(x_k_minus(2,1),x_k_minus(1,1));hypot(x_k_minus(3,1),x_k_minus(4,1));atan2d(x_k_minus(4,1),x_k_minus(3,1))];
    
    %Calculate measurement mapping Jacobian 
    range_denom=1/hypot(x_k_minus(1,1),x_k_minus(2,1));
    range_dx=x_k_minus(1,1)*range_denom;
    range_dy=x_k_minus(2,1)*range_denom;
      
    %With position
    azi_dx=-1*x_k_minus(2,1)/(x_k_minus(1,1)^2+x_k_minus(2,1)^2);
    azi_dy=1/(x_k_minus(1,1)+(x_k_minus(2,1)^2/x_k_minus(1,1)));
    
    stw_denom=1/hypot(x_k_minus(3,1),x_k_minus(4,1));
    stw_dx=x_k_minus(3,1)*stw_denom;
    stw_dy=x_k_minus(4,1)*stw_denom;
    
    hdg_dx=-1*x_k_minus(4,1)/(x_k_minus(3,1)^2+x_k_minus(4,1)^2);
    hdg_dy=1/(x_k_minus(3,1)+(x_k_minus(4,1)^2/x_k_minus(3,1)));
%     speed_denom=1/hypot((x_k_minus(3,1)-x_k_minus(5,1)),(x_k_minus(4,1)-x_k_minus(6,1)));
%     speed_dx=(x_k_minus(3,1)-x_k_minus(5,1))*speed_denom;
%     speed_dy=(x_k_minus(4,1)-x_k_minus(6,1))*speed_denom;
%     
%     hdg_dx=-1*(x_k_minus(4,1)-x_k_minus(6,1))/((x_k_minus(3,1)-x_k_minus(5,1))^2+(x_k_minus(4,1)-x_k_minus(6,1))^2);
%     hdg_dy=1/((x_k_minus(3,1)-x_k_minus(5,1))+((x_k_minus(4,1)-x_k_minus(6,1))^2/(x_k_minus(3,1)-x_k_minus(5,1))));
    
    H=[range_dx range_dy 0 0 0 0;azi_dx azi_dy 0 0 0 0;0 0 stw_dx stw_dy 0 0;0 0 hdg_dx hdg_dy 0 0];
    
    
    %Calculate innovation
    r_k=z_k-h;
       
    %Calculate innovation covariance
    S=H*P_minus*H'+R;
    
    %Update Kalman Gain after estimate vector calculation
    K=(P_minus*H')/S;
    
    
    %Update state vector
    x_k_plus=x_k_minus+K*r_k;
    
    
        
    %Update error covariance matrix with new Kalman Gain
    %P_plus=(eye(4)-K*H)*P_minus;
    
    %Joseph stabilized equation
    P_plus=(eye(6)-K*H)*P_minus*(eye(6)-K*H)'+K*R*K';
    else
        x_k_plus=x_k_minus;
        P_plus=P_minus;
    end
    x_err=x_act-x_k_plus;
    epsilon=x_err.'*inv(P_plus)*x_err;
    data(7,ii)=x_k_plus(1,1);
    data(8,ii)=x_k_plus(2,1);
    data(9,ii)=x_k_plus(3,1);
    data(10,ii)=x_k_plus(4,1);
    data(11,ii)=x_k_plus(5,1);
    data(12,ii)=x_k_plus(6,1);
    data(13,ii)=hypot(x_k_plus(3,1),x_k_plus(4,1));
    data(14,ii)=atand(x_k_plus(4,1)/x_k_plus(3,1));
    data(15,ii)=range_x;
    data(16,ii)=range_y;
    data(17,ii)=norm(P_plus);
    data(18,ii)=range-hypot(x_k_plus(1,1),x_k_plus(2,1));
    data(19,ii)=azi-atan2d(x_k_plus(2,1),x_k_plus(1,1));
    data(20,ii)=S(1,1);
    data(21,ii)=S(2,1);
    data(22,ii)=S(3,1);
    data(23,ii)=S(4,1);
    data(24,ii)=r_k(1,1);
    data(25,ii)=r_k(2,1);
    data(26,ii)=r_k(3,1);
    data(27,ii)=r_k(4,1);
    data(28,ii)=epsilon;
end


% avg(1,1)=hypot(mean(data(11,1600*freq:2000*freq)),mean(data(12,1600*freq:2000*freq)));
% avg(2,1)=atan2d(mean(data(12,1600*freq:2000*freq)),mean(data(11,1600*freq:2000*freq)));
% avg(3,1)=hypot(mean(data(11,3600*freq:steps)),mean(data(12,3600*freq:steps)));
% avg(4,1)=atan2d(mean(data(12,3600*freq:steps)),mean(data(11,3600*freq:steps)));
mag=hypot(data(11,:),data(12,:));
dir=atan2d(data(12,:),data(11,:));

figure
plot(1:steps,data(17,1:steps))
title('Error Covariance Matrix Norm vs Filter Step')
xlabel('Filter Step')
ylabel('Error Covariance Matrix Norm')

figure
plot(1:steps,data(18,:))
%axis([0 steps -5 5])
title('Filter Error vs Filter Step')
xlabel('Filter Step')
ylabel('Filter Error (m)')

figure
plot(1:steps,data(19,:))
%axis([0 steps -5 5])
title('Filter Error vs Filter Step')
xlabel('Filter Step')
ylabel('Filter Error (degrees)')

figure
plot(data(7,:),data(8,:),data(15,:),data(16,:))
axis equal
title('Acutal Position vs Filter Output Position')
xlabel('x position (m)')
ylabel('y position (m)')
legend('Filter Output','Actual Position','location','southeast')

figure
yyaxis left
plot(1:steps,data(24,:),1:steps,data(25,:))
yyaxis right
plot(1:steps,data(26,:),1:steps,data(27,:))

figure
plot(1:steps,set,1:steps,mag)

figure
plot(1:steps,drift,1:steps,dir)

figure
plot(data(28,:))

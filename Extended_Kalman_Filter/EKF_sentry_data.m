clear variables;
load D:\Documents\MATLAB\Research_Project\lab_soln\kf.mat
%Define # of filter steps desired
%steps=3300;

%Define time step size
dt=1;

%Initialize data matrix
data=zeros(24,length(kf.t)-1);
avg=zeros(6,1);
%rts_out=zeros(6,steps-5);

%Define boundary conditions
x0_pos=kf.dr_ahrs.e(1);
y0_pos=kf.dr_ahrs.n(1);
%stw=.50;
%hdg=30;
%set=.3;
%drift=30;
x0_vel=0;
y0_vel=0;

%Initialize estimated state vector
%x_k_minus=[x0_pos;y0_pos;x0_vel;y0_vel;set*cosd(drift);set*sind(drift)];
x_k_minus=[x0_pos;y0_pos;x0_vel;y0_vel;0;0];
%Compute Q, the covariance matrix for noise associated with the state vector
var_Q=.3^2;
%var_Q=.001^2;
%G=[.5*(dt^2); .5*(dt^2);dt; dt;dt;dt];
%Q=G*var_Q*G';
Q=diag([var_Q,var_Q,var_Q,var_Q,var_Q,var_Q]);

%Compute C, the linear transformation matrix for white noise
%[U_P,T_P,V_P]=svd(Q);
%S_P=sqrt(T_P);
%C=U_P*S_P;

%Q with different sigma values
%var_Qp=.001^2;
%var_Qs=.001^2;
%Q=diag([var_Qp,var_Qp,var_Qs,var_Qs,var_Qs,var_Qs]);

%Compute R, the covariance matrix associated with measurement error
var_Rr=10^2;
var_Ra=1^2;
%var_Rr=.001^2;
%var_Ra=.001^2;
%R=diag([var_Rr,var_Ra]);
R=var_Rr;

%Initialize state transition matrix
F=[1 0 0 0 dt 0;0 1 0 0 0 dt;0 0 0 0 1 0;0 0 0 0 0 1;0 0 0 0 1 0;0 0 0 0 0 1];

%Initialize error covariance matrix
var_Pr=1^2;
var_Ps=.2^2;
P_minus=diag([var_Pr,var_Pr,var_Ps,var_Ps,var_Ps,var_Ps]);
P_plus=P_minus;

%Initialize measurement vector
z_k=kf.rng.hr(1);

%Initialize state vector
x_k_plus=x_k_minus;

%jj=1;
for ii=2:length(kf.t)
        
    %Calculate state transition matrix
    dt=kf.t(ii)-kf.t(ii-1);
    F=[1 0 0 0 dt 0;0 1 0 0 0 dt;0 0 0 0 1 0;0 0 0 0 0 1;0 0 0 0 1 0;0 0 0 0 0 1];
    
    %Calculate estimate for next step
    u_k=[kf.dr_ahrs.de(ii);kf.dr_ahrs.dn(ii);kf.dr_ahrs.de(ii)/dt;kf.dr_ahrs.dn(ii)/dt;0;0];
    x_k_minus=F*x_k_plus+u_k;
    
    %System projection with linear constraint (cog)
    %if ii>3
    %    dx=data(7,ii-1)-data(7,ii-3);
    %    dy=data(8,ii-1)-data(8,ii-3);
    %    cog=atan2d(dy,dx);
    %    D=[0 0 -1*tand(cog) 1 0 0];
    %    d=0;
    %    [U_Q,S_Q,V_Q]=svd(D');
    %    U_Q(:,1)=[];
    %    N=U_Q*U_Q';
    %    Q=N*Q*N;
    %end
    
    %System projection with non-linear constraint (stw)
    %g_x=hypot((x_k_minus(3,1)-x_k_minus(5,1)),(x_k_minus(4,1)-x_k_minus(6,1)));
    %dg_dx=(x_k_minus(3,1)-x_k_minus(5,1))/g_x;
    %dg_dy=(x_k_minus(4,1)-x_k_minus(6,1))/g_x;
    %D=[0 0 dg_dx dg_dy -1*dg_dx -1*dg_dy];
    %d=hypot(u_k(3,1),u_k(4,1))-g_x+(D*x_k_minus);
    %[U_Q,S_Q,V_Q]=svd(D');
    %U_Q(:,1)=[];
    %N=U_Q*U_Q';
    %Q=N*Q*N;
    
    %Calculate error covariance matrix for next step
    P_minus=F*P_plus*F'+Q;
        
    data(1,ii-1)=x_k_minus(1,1);
    data(2,ii-1)=x_k_minus(2,1);
    data(3,ii-1)=x_k_minus(3,1);
    data(4,ii-1)=x_k_minus(4,1);
    data(5,ii-1)=x_k_minus(5,1);
    data(6,ii-1)=x_k_minus(6,1);
    
    %Calculate measurement vector
    z_k=kf.rng.hr(ii);
    
    %Calculate non-linear measurement
    %h_ang=atan2d(kf.rng.n(ii)-x_k_minus(2,1),kf.rng.e(ii)-x_k_minus(1,1));
    h=hypot(kf.rng.e(ii)-x_k_minus(1,1),kf.rng.n(ii)-x_k_minus(2,1));
    %h=sqrt((kf.rng.e(ii)-x_k_minus(1,1))^2+(kf.rng.n(ii)-x_k_minus(2,1))^2);
   
    %Calculate measurement mapping Jacobian 
    range_denom=1/hypot(kf.rng.e(ii)-x_k_minus(1,1),kf.rng.n(ii)-x_k_minus(2,1));
    range_dx=kf.rng.e(ii)-x_k_minus(1,1)*range_denom;
    range_dy=kf.rng.n(ii)-x_k_minus(2,1)*range_denom;
   
    %With position
    %azi_dx=-1*kf.rng.n(ii)-x_k_minus(2,1)/((kf.rng.e(ii)-x_k_minus(1,1))^2+(kf.rng.n(ii)-x_k_minus(2,1))^2);
    %azi_dy=1/((kf.rng.e(ii)-x_k_minus(1,1))+((kf.rng.n(ii)-x_k_minus(2,1))^2/(kf.rng.e(ii)-x_k_minus(1,1))));
    H=[range_dx range_dy 0 0 0 0];
    
    %Calculate innovation
    r_k=z_k-h;
       
    %Calculate innovation covariance
    S=H*P_minus*H'+R;
    
    %Update Kalman Gain after estimate vector calculation
    K=(P_minus*H')/S;
    
    %Gain projection with linear constraint (cog)
    %if ii>3
    %    dx=data(7,ii-1)-data(7,ii-3);
    %    dy=data(8,ii-1)-data(8,ii-3);
    %    cog=atan2d(dy,dx);
    %    D=[0 0 -1*tand(cog) 1 0 0];
    %    d=0;
    %    K=K-D'/(D*D')*(D*x_k_minus-d)/(r_k'/S*r_k)*(r_k'/S);
    %end
    
    %Update state vector
    x_k_plus=x_k_minus+K*r_k;
        
    %Update error covariance matrix with new Kalman Gain
    %P_plus=(eye(4)-K*H)*P_minus;
    
    %Joseph stabilized equation
    P_plus=(eye(6)-K*H)*P_minus*(eye(6)-K*H)'+K*R*K';
    
    data(7,ii-1)=x_k_plus(1,1);
    data(8,ii-1)=x_k_plus(2,1);
    data(9,ii-1)=x_k_plus(3,1);
    data(10,ii-1)=x_k_plus(4,1);
    data(11,ii-1)=x_k_plus(5,1);
    data(12,ii-1)=x_k_plus(6,1);
    data(13,ii-1)=z_k;
    data(14,ii-1)=h;
    %data(13,ii-1)=hypot(x_k_plus(3,1),x_k_plus(4,1));
    %data(14,ii-1)=atand(x_k_plus(4,1)/x_k_plus(3,1));
    %data(15,ii-1)=range_x;
    %data(16,ii-1)=range_y;
    data(17,ii-1)=norm(P_plus);
    %data(18,ii-1)=range-hypot(x_k_plus(1,1),x_k_plus(2,1));
    %data(19,ii-1)=azi-atan2d(x_k_plus(2,1),x_k_plus(1,1));
    data(20,ii-1)=S(1,1);
    %data(21,ii-1)=S(2,1);
    data(22,ii-1)=r_k(1,1);
    %data(23,ii-1)=r_k(2,1);
end

%avg(1,1)=mean(data(11,1:2000));
%avg(2,1)=mean(data(12,1:2000));
%avg(3,1)=atan2d(mean(data(12,1:2000)),mean(data(11,1:2000)));
%avg(4,1)=mean(data(11,2001:steps));
%avg(5,1)=mean(data(12,2001:steps));
%avg(6,1)=atan2d(mean(data(12,2001:steps)),mean(data(11,2001:steps)));

figure
plot(1:length(kf.t)-1,data(17,1:length(kf.t)-1))
title('Error Covariance Matrix Norm vs Filter Step')
xlabel('Filter Step')
ylabel('Error Covariance Matrix Norm')

%figure
%plot(1:length(kf.t)-1),data(18,:))
%axis([0 steps -5 5])
%title('Filter Error vs Filter Step')
%xlabel('Filter Step')
%ylabel('Filter Error (m)')

figure
plot(data(7,:),data(8,:))%,data(15,:),data(16,:))
%title('Acutal Position vs Filter Output Position')
%xlabel('x position (m)')
%ylabel('y position (m)')
%legend('Filter Output','Actual Position','location','southeast')

%figure
%plot(rts_out(1,:),rts_out(2,:),data(15,6:steps),data(16,6:steps))
%title('Acutal Position vs RTS Output Position')
%xlabel('x position (m)')
%ylabel('y position (m)')
%legend('RTS Output','Actual Position','location','southeast')


    

clear variables;

%Define time step size
freq=10;
dt=1/freq;

%Define # of filter steps desired
steps=4000*freq;

%Initialize data matrix
data=zeros(24,steps);
avg=zeros(4,1);

%Define boundary conditions
x0_pos=0;
y0_pos=0;
stw=.50;
hdg=30;
set=.3;
drift=60;
x0_vel=stw*cosd(hdg);%+set*cosd(drift);
y0_vel=stw*sind(hdg);%+set*sind(drift);
S=[0;0];
r_k=[0;0];

%Initialize estimated state vector
%x_k_minus=[x0_pos;y0_pos;x0_vel;y0_vel;set*cosd(drift);set*sind(drift)];
x_k_minus=[x0_pos;y0_pos;x0_vel;y0_vel;0;0];
%Compute Q, the covariance matrix for noise associated with the state vector
var_Q=diag(.3^2,.3^2,.3^2,.3^2,.025^2,.025^2);
%var_Q=.0001^2;
G=[.5*(dt^2); .5*(dt^2);dt; dt;dt;dt];
Q=G*G'*var_Q;

%Compute C, the linear transformation matrix for white noise
[U_P,T_P,V_P]=svd(Q);
S_P=sqrt(T_P);
C=U_P*S_P;

%Compute R, the covariance matrix associated with measurement error
var_Rr=8^2;
var_Ra=1^2;
%var_Rr=.0001^2;
%var_Ra=.0001^2;
R=diag([var_Rr,var_Ra]);

%Initialize error covariance matrix
var_Pr=10^2;
var_Ps=1^2;
var_Pa=1^2;
P_minus=diag([var_Pr,var_Pr,var_Ps,var_Ps,var_Ps,var_Ps]);
P_plus=P_minus;

%Initialize measurement vector
range_x=0;
range_y=0;
range=0;
azi=0;
z_k=[0;0];

%Initialize state vector
x_k_plus=x_k_minus;

%F=[1 0 0 0 dt 0 0 0;0 1 0 0 0 dt 0 0;0 0 0 0 1 0 0 0;0 0 0 0 0 1 0 0;0 0 0 0 1 0 0 0;0 0 0 0 0 1 0 0;0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0];

for ii=1:steps
    %Define courses for simulation
    if ii<=400*freq
        stw=.5;
        hdg=30;
    elseif ii<=800*freq
            stw=.5;
            hdg=120;
    elseif ii<=1200*freq
            stw=.5;
            hdg=210;
    elseif ii<=1600*freq
            stw=.5;
            hdg=120;
    elseif ii<=2400*freq
        stw=.5;
        hdg=30;
    elseif ii<=2800*freq
        stw=.5;
        hdg=120;
    elseif ii<=3200*freq
        stw=.5;
        hdg=210;
    elseif ii<=3600*freq
        stw=.5;
        hdg=120;
    else
        stw=.5;
        hdg=30;
    end
    
    %Specify currents for simulation
    if ii<=2000*freq
        set=.3;
        drift=60;
    else
        set=.3;
        drift=30;
    end
    
    if ii==2000*freq+1
        x_k_plus(3,1)=stw*cosd(hdg);
        x_k_plus(4,1)=stw*sind(hdg);
        x_k_plus(5,1)=0;
        x_k_plus(6,1)=0;
        %P_plus(3,3)=.2^2;
        %P_plus(4,4)=.2^2;
        %P_plus(5,5)=.2^2;
        %P_plus(6,6)=.2^2;
        
        %P_plus=diag([var_Pr,var_Pr,var_Ps,var_Ps,var_Ps,var_Ps,var_Pa,var_Pa]);
        %P_plus=P_plus.*15;
        
        P_plus(3,3)=P_plus(3,3).*1.5;
        P_plus(4,4)=P_plus(4,4).*1.5;
        P_plus(5,5)=P_plus(5,5).*1.5;
        P_plus(6,6)=P_plus(6,6).*1.5;
    end
    x_vel=stw*cosd(hdg)+set*cosd(drift);
    y_vel=stw*sind(hdg)+set*sind(drift);
    %sog=hypot(x_vel,y_vel);
    %cog=atan2d(y_vel,x_vel);
    
    
    %Calculate noisy estimate vector for next step
    u_k=[stw*dt*cosd(hdg);stw*dt*sind(hdg);stw*cosd(hdg);stw*sind(hdg);0;0;0;0];
    w_k=randn(8,1);
    w_k=C*w_k;
    
    
    
    f=[x_k_plus(1,1)+(x_k_plus(5,1)*dt); x_k_plus(2,1)+(x_k_plus(6,1)*dt); new_x; new_y; new_x; new_y];
    x_k_minus=f+u_k+w_k;
    
    %x_k_minus=F*x_k_plus+u_k+w_k;
   
    
   
    
    
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
    
    
    if mod(ii,freq)==0
    v_k=randn(2,1);
    [V_R,D_R]=eig(R);
    [~,ind] = sort(diag(D_R),'descend');
    D_Rs = D_R(ind,ind);
    V_Rs = V_R(:,ind);
    v_k=V_Rs*D_Rs*v_k;
    z_k=[range;azi]+v_k;
    
    %Calculate non-linear measurement
    h_ang=atan2d(x_k_minus(2,1),x_k_minus(1,1));
    h=[hypot(x_k_minus(1,1),x_k_minus(2,1));h_ang];
    
    %Calculate measurement mapping Jacobian 
    range_denom=1/hypot(x_k_minus(1,1),x_k_minus(2,1));
    range_dx=x_k_minus(1,1)*range_denom;
    range_dy=x_k_minus(2,1)*range_denom;
      
    %With position
    azi_dx=-1*x_k_minus(2,1)/(x_k_minus(1,1)^2+x_k_minus(2,1)^2);
    azi_dy=1/(x_k_minus(1,1)+(x_k_minus(2,1)^2/x_k_minus(1,1)));
    H=[range_dx range_dy 0 0 0 0 0 0;azi_dx azi_dy 0 0 0 0 0 0];
    
    
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
    P_plus=(eye(8)-K*H)*P_minus*(eye(8)-K*H)'+K*R*K';
    else
        x_k_plus=x_k_minus;
        P_plus=P_minus;
    end
    
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
    data(22,ii)=r_k(1,1);
    data(23,ii)=r_k(2,1);
end


avg(1,1)=hypot(mean(data(11,1600*freq:2000*freq)),mean(data(12,1600*freq:2000*freq)));
avg(2,1)=atan2d(mean(data(12,1600*freq:2000*freq)),mean(data(11,1600*freq:2000*freq)));
avg(3,1)=hypot(mean(data(11,3600*freq:steps)),mean(data(12,3600*freq:steps)));
avg(4,1)=atan2d(mean(data(12,3600*freq:steps)),mean(data(11,3600*freq:steps)));


figure
plot(1:steps,data(17,1:steps))
title('Error Covariance Matrix Norm vs Filter Step')
xlabel('Filter Step')
ylabel('Error Covariance Matrix Norm')

figure
plot(1:steps,data(18,:))
axis([0 steps -5 5])
title('Filter Error vs Filter Step')
xlabel('Filter Step')
ylabel('Filter Error (m)')

figure
plot(1:steps,data(19,:))
axis([0 steps -5 5])
title('Filter Error vs Filter Step')
xlabel('Filter Step')
ylabel('Filter Error (degrees)')

figure
plot(data(7,:),data(8,:),data(15,:),data(16,:))
title('Acutal Position vs Filter Output Position')
xlabel('x position (m)')
ylabel('y position (m)')
legend('Filter Output','Actual Position','location','southeast')


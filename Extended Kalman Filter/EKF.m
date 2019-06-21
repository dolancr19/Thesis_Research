clear variables;

%Define # of filter steps desired
steps=10800;

%Define time step size
dt=1;

%Initialize data matrix
data=zeros(23,steps);

avg=zeros(6,1);

%Define boundary conditions
x0_pos=0;
y0_pos=0;
stw=.50;
hdg=30;
set=.2;
drift=30;
x0_vel=stw*cosd(hdg)+set*cosd(drift);
y0_vel=stw*sind(hdg)+set*sind(drift);
sog=hypot(x0_vel,y0_vel);
cog=atan2d(y0_vel,x0_vel);

%Initialize estimated state vector
x_k_minus=[x0_pos;y0_pos;x0_vel;y0_vel;set*cosd(drift);set*sind(drift)];

%Compute Q, the covariance matrix for noise associated with the state vector
var_Q=.02^2;
G=[.5*(dt^2); .5*(dt^2);dt; dt;dt;dt];
Q=G*var_Q*G';

%Compute C, the linear transformation matrix for white noise
[U,T,V]=svd(Q);
S=sqrt(T);
C=U*S;

%Q with different sigma values
%var_Qp=.001^2;
%var_Qs=.001^2;
%Q=diag([var_Qp,var_Qp,var_Qs,var_Qs,var_Qs,var_Qs]);

%Compute R, the covariance matrix associated with measurement error
var_Rr=.02^2;
var_Ra=2^2;
R=diag([var_Rr,var_Ra]);

%Initialize state transition matrix
F=[1 0 0 0 dt 0;0 1 0 0 0 dt;0 0 0 0 1 0;0 0 0 0 0 1;0 0 0 0 1 0;0 0 0 0 0 1];

%Initialize error covariance matrix
var_P=10^2;
P_minus=diag([var_P,var_P,var_P,var_P,var_P,var_P]);

P_plus=P_minus;

%Initialize measurement vector
range_x=0;
range_y=0;
range=0;
azi=0;

z_k=[0;0];

%Initialize state vector
x_k_plus=x_k_minus;

jj=1;
for ii=1:steps
    
    if ii<=3600
        stw=.5;
        hdg=30;
    elseif ii<=7200
            stw=.5;
            hdg=90;
    else
        stw=.5;
        hdg=210;        
    end
    
    x_vel=stw*cosd(hdg)+set*cosd(drift);
    y_vel=stw*sind(hdg)+set*sind(drift);
    sog=hypot(x_vel,y_vel);
    cog=atan2d(y_vel,x_vel);
    
    %Calculate noisy estimate vector for next step
    u_k=[stw*dt*cosd(hdg);stw*dt*sind(hdg);stw*cosd(hdg);stw*sind(hdg);0;0];
    
    w_k=randn(6,1);
    %[V_Q,D_Q]=eig(Q);
    %[dq,indq] = sort(diag(D_Q),'descend');
    %D_Qs = D_Q(indq,indq);
    %V_Qs = V_Q(:,indq);
    
    %w_k=V_Qs*D_Qs*w_k;
    w_k=C*w_k;
    x_k_minus=F*x_k_plus+u_k+w_k;
    
    %Calculate error covariance matrix for next step
    P_minus=F*P_plus*F'+Q;
    
    
    data(1,ii)=x_k_minus(1,1);
    data(2,ii)=x_k_minus(2,1);
    data(3,ii)=x_k_minus(3,1);
    data(4,ii)=x_k_minus(4,1);
    data(5,ii)=x_k_minus(5,1);
    data(6,ii)=x_k_minus(6,1);
    
    %Calculate measurement mapping Jacobian 
    range_denom=1/hypot(x_k_minus(1,1),x_k_minus(2,1));
    range_dx=x_k_minus(1,1)*range_denom;
    range_dy=x_k_minus(2,1)*range_denom;
   
    %With velocity
    %cog_dx=-1*x_k_minus(4,1)/(x_k_minus(3,1)^2+x_k_minus(4,1)^2);
    %cog_dy=1/(x_k_minus(3,1)+(x_k_minus(4,1)^2/x_k_minus(3,1)));
    %H=[range_dx range_dy 0 0 0 0;0 0 cog_dx cog_dy 0 0];
    
    %With position
    cog_dx=-1*x_k_minus(2,1)/(x_k_minus(1,1)^2+x_k_minus(2,1)^2);
    cog_dy=1/(x_k_minus(1,1)+(x_k_minus(2,1)^2/x_k_minus(1,1)));
    H=[range_dx range_dy 0 0 0 0;cog_dx cog_dy 0 0 0 0];
    
    
    
    %Update Kalman Gain after estimate vector calculation
    K=(P_minus*H')/(H*P_minus*H'+R);
    s=H*P_minus*H'+R;
    %Calculate noisy measurement vector
    range_x=range_x+x_vel*dt;
    range_y=range_y+y_vel*dt;
    range=hypot(range_x,range_y);
    azi=atan2d(range_y,range_x);
    
    %sigma_R=.02*range;
    %if sigma_R<7
    %    var_Rr=sigma_R^2;
    %else
    %    var_Rr=7^2;
    %end
    var_Rr=7^2;
    var_Ra=2^2;
    R=diag([var_Rr,var_Ra]);
    v_k=randn(2,1);
    [V_R,D_R]=eig(R);
    [d,ind] = sort(diag(D_R),'descend');
    D_Rs = D_R(ind,ind);
    V_Rs = V_R(:,ind);
    
    v_k=V_Rs*D_Rs*v_k;
    z_k=[range;azi]+v_k;
            
    %Update state vector
    
    %With velocity
    %if x_k_minus(4,1)>0
    %    if x_k_minus(3,1)>0
    %        h_ang=atand(x_k_minus(4,1)/x_k_minus(3,1));
    %    else
    %        h_ang=180+atand(x_k_minus(4,1)/x_k_minus(3,1));
    %    end
    %else
    %    if x_k_minus(3,1)>0
    %        h_ang=atand(x_k_minus(4,1)/x_k_minus(3,1));
    %    else
    %        h_ang=180+atand(x_k_minus(4,1)/x_k_minus(3,1));
    %    end
    %end
    
    %With position
    %if x_k_minus(2,1)>0
    %    if x_k_minus(1,1)>0
    %        h_ang=atand(x_k_minus(2,1)/x_k_minus(1,1));
    %    else
    %        h_ang=180+atand(x_k_minus(2,1)/x_k_minus(1,1));
    %    end
    %else
    %    if x_k_minus(1,1)>0
    %        h_ang=atand(x_k_minus(2,1)/x_k_minus(1,1));
    %    else
    %        h_ang=180+atand(x_k_minus(2,1)/x_k_minus(1,1));
    %    end
    %end
    
    %With atan2d
    h_ang=atan2d(x_k_minus(2,1),x_k_minus(1,1));
    
    h=[hypot(x_k_minus(1,1),x_k_minus(2,1));h_ang];
    x_k_plus=x_k_minus+K*(z_k-h);
    nu=z_k-h;
    %if jj<= 50
    %    remain = mod(ii,60);
    %    if remain == 0
    %        x_k_plus(1,1)=ii*sog*dt*cosd(cog);
    %        x_k_plus(2,1)=ii*sog*dt*sind(cog);
    %        jj=jj+1;
    %    end
    %end
    
            
    
    %Update error covariance matrix with new Kalman Gain
    %P_plus=(eye(4)-K*H)*P_minus;
    
    %Joseph stabilized equation
    P_plus=(eye(6)-K*H)*P_minus*(eye(6)-K*H)'+K*R*K';
    
    remain=mod(ii,300);
    %if remain==0
    %    x_k_plus(3,1)=mean(data(9,:));
    %    x_k_plus(4,1)=mean(data(10,:));
    %    x_k_plus(5,1)=mean(data(11,:));
     %   x_k_plus(6,1)=mean(data(12,:));
    %end
    
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
    data(19,ii)=cog-atan2d(x_k_plus(2,1),x_k_plus(1,1));
    data(20,ii)=s(1,1);
    data(21,ii)=s(2,1);
    data(22,ii)=nu(1,1);
    data(23,ii)=nu(2,1);
    
    
    
end
avg(1,1)=mean(data(9,:));
avg(2,1)=mean(data(10,:));
avg(3,1)=atan2d(mean(data(10,:)),mean(data(9,:)));
avg(4,1)=mean(data(11,:));
avg(5,1)=mean(data(12,:));
avg(6,1)=atan2d(mean(data(12,:)),mean(data(11,:)));

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
plot(data(7,:),data(8,:),data(15,:),data(16,:))
title('Acutal Position vs Filter Output Position')
xlabel('x position (m)')
ylabel('y position (m)')
legend('Filter Output','Actual Position','location','southeast')

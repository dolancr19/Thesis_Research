clear variables;

%Define # of filter steps desired
steps=3300;

%Define time step size
dt=1;

%Initialize data matrix
data=zeros(24,steps);
avg=zeros(6,1);
rts_out=zeros(6,steps-5);

%Define boundary conditions
x0_pos=0;
y0_pos=0;
stw=.50;
hdg=30;
set=.3;
drift=30;
x0_vel=stw*cosd(hdg);%+set*cosd(drift);
y0_vel=stw*sind(hdg);%+set*sind(drift);

%Initialize estimated state vector
%x_k_minus=[x0_pos;y0_pos;x0_vel;y0_vel;set*cosd(drift);set*sind(drift)];
x_k_minus=[x0_pos;y0_pos;x0_vel;y0_vel;0;0];
%Compute Q, the covariance matrix for noise associated with the state vector
var_Q=.015^2;
%var_Q=.001^2;
G=[.5*(dt^2); .5*(dt^2);dt; dt;dt;dt];
Q=G*var_Q*G';

%Compute C, the linear transformation matrix for white noise
[U_P,T_P,V_P]=svd(Q);
S_P=sqrt(T_P);
C=U_P*S_P;

%Q with different sigma values
%var_Qp=.001^2;
%var_Qs=.001^2;
%Q=diag([var_Qp,var_Qp,var_Qs,var_Qs,var_Qs,var_Qs]);

%Compute R, the covariance matrix associated with measurement error
var_Rr=5^2;
var_Ra=1^2;
%var_Rr=.001^2;
%var_Ra=.001^2;
R=diag([var_Rr,var_Ra]);

%Initialize state transition matrix
F=[1 0 0 0 dt 0;0 1 0 0 0 dt;0 0 0 0 1 0;0 0 0 0 0 1;0 0 0 0 1 0;0 0 0 0 0 1];

%Initialize error covariance matrix
var_Pr=1^2;
var_Ps=.2^2;
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

%jj=1;
for ii=1:steps
    %Define courses for simulation
    if ii<=600
        stw=.5;
        hdg=30;
    elseif ii<=900
            stw=.5;
            hdg=120;
    elseif ii<=1500
            stw=.5;
            hdg=210;
    elseif ii<=1800
            stw=.5;
            hdg=120;
    elseif ii<=2400
        stw=.5;
        hdg=30;
    elseif ii<=2700
        stw=.5;
        hdg=120;
    else
        stw=.5;
        hdg=210;
    end
    
    %Specify currents for simulation
    if ii<=2000
        set=.3;
        drift=30;
    else
        set=.5;
        drift=30;
    end
    
    %if ii==2001
    %    f_k_plus(3,1)=stw*cosd(hdg);
    %    f_k_plus(4,1)=stw*sind(hdg);
    %    f_k_plus(5,1)=0;
    %    f_k_plus(6,1)=0;
        %P_plus(3,3)=.2^2;
        %P_plus(4,4)=.2^2;
        %P_plus(5,5)=.2^2;
        %P_plus(6,6)=.2^2;
        
    %    P_plus=diag([P_plus(1,1),P_plus(2,2),.2^2,.2^2,.2^2,.2^2]);
        %P_plus=P_plus.*10;
    %end
    x_vel=stw*cosd(hdg)+set*cosd(drift);
    y_vel=stw*sind(hdg)+set*sind(drift);
    %sog=hypot(x_vel,y_vel);
    %cog=atan2d(y_vel,x_vel);
    
    
    %Calculate noisy estimate vector for next step
    u_k=[stw*dt*cosd(hdg);stw*dt*sind(hdg);stw*cosd(hdg);stw*sind(hdg);0;0];
    w_k=randn(6,1);
    w_k=C*w_k;
    
    %With different sigma values for each measurement
    %[V_Q,D_Q]=eig(Q);
    %[dq,indq] = sort(diag(D_Q),'descend');
    %D_Qs = D_Q(indq,indq);
    %V_Qs = V_Q(:,indq);
    %w_k=V_Qs*D_Qs*w_k;
    
    x_k_minus=F*x_k_plus+u_k+w_k;
    
    %if ii<=6
    %    field="xm"+ii;
    %    rts.(field)=x_k_minus;
    %else
    %    rts.xm1=rts.xm2;
    %    rts.xm2=rts.xm3;
    %    rts.xm3=rts.xm4;
    %    rts.xm4=rts.xm5;
    %    rts.xm5=rts.xm6;
    %    rts.xm6=x_k_minus;
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
    
    %System projection with linear constraint (cog)
    %Q=diag([var_Q,var_Q,var_Q,var_Q,var_Q,var_Q]);
    if ii>2
        dx=data(7,ii-1)-data(7,ii-2);
        dy=data(8,ii-1)-data(8,ii-2);
        cog=atan2d(dy,dx);
        D=[0 0 -1*tand(cog) 1 0 0];
        d=0;
        [U_Q,S_Q,V_Q]=svd(D');
        U_Q(:,1)=[];
        N=U_Q*U_Q';
        Q=N*Q*N;
    end
    
    
    %Calculate error covariance matrix for next step
    P_minus=F*P_plus*F'+Q;
    %if ii<=6
    %    field="pm"+ii;
    %    rts.(field)=P_minus;
    %else
    %    rts.pm1=rts.pm2;
    %    rts.pm2=rts.pm3;
    %    rts.pm3=rts.pm4;
    %    rts.pm4=rts.pm5;
    %    rts.pm5=rts.pm6;
    %    rts.pm6=P_minus;
    %end
        
    
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
    
    %Determine error in range based on the range value
    %Tuning will be more complicated than the simple equation below.
    %Shelved for now.
    %sigma_R=.02*range;
    %if sigma_R<7
    %    var_Rr=sigma_R^2;
    %else
    %    var_Rr=7^2;
    %end
    %R=diag([var_Rr,var_Ra]);
    
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
   
    %With velocity
    %Only works with a single course, not useful
    %cog_dx=-1*x_k_minus(4,1)/(x_k_minus(3,1)^2+x_k_minus(4,1)^2);
    %cog_dy=1/(x_k_minus(3,1)+(x_k_minus(4,1)^2/x_k_minus(3,1)));
    %H=[range_dx range_dy 0 0 0 0;0 0 cog_dx cog_dy 0 0];
    
    %With position
    azi_dx=-1*x_k_minus(2,1)/(x_k_minus(1,1)^2+x_k_minus(2,1)^2);
    azi_dy=1/(x_k_minus(1,1)+(x_k_minus(2,1)^2/x_k_minus(1,1)));
    H=[range_dx range_dy 0 0 0 0;azi_dx azi_dy 0 0 0 0];
    
    %Calculate innovation
    r_k=z_k-h;
       
    %Calculate innovation covariance
    S=H*P_minus*H'+R;
    
    %Update Kalman Gain after estimate vector calculation
    K=(P_minus*H')/S;
    
    %Gain projection with non-linear constraint (stw)
    %g_x=hypot((x_k_minus(3,1)-x_k_minus(5,1)),(x_k_minus(4,1)-x_k_minus(6,1)));
    %dg_dx=(x_k_minus(3,1)-x_k_minus(5,1))/g_x;
    %dg_dy=(x_k_minus(4,1)-x_k_minus(6,1))/g_x;
    %D=[0 0 dg_dx dg_dy -1*dg_dx -1*dg_dy];
    %d=hypot(u_k(3,1),u_k(4,1))-g_x+(D*x_k_minus);
    %K=K-D'/(D*D')*(D*x_k_minus-d)/(r_k'/S*r_k)*(r_k'/S);   
    
    %Gain projection with non-linear constraint (range)
    %g_x=x_k_minus(1,1)/cosd(z_k(2,1));
    %dg_dx=1/cosd(z_k(2,1));
    %D=[dg_dx 0 0 0 0 0];
    %d=z_k(1,1)-g_x+(D*x_k_minus);
    %K=K-D'/(D*D')*(D*x_k_minus-d)/(r_k'/S*r_k)*(r_k'/S);
    
    %Gain projection with non-linear constraint (azi)
    %g_x=x_k_minus(2,1)/x_k_minus(1,1);
    %dg_dx=-1*x_k_minus(2,1)/(x_k_minus(1,1)^2);
    %dg_dy=1/x_k_minus(1,1);
    %D=[dg_dx dg_dy 0 0 0 0];
    %d=tand(z_k(2,1))-g_x+(D*x_k_minus);
    %K=K-D'/(D*D')*(D*x_k_minus-d)/(r_k'/S*r_k)*(r_k'/S);
    
    %Gain projection with linear constraint (cog)
    if ii>2
        dx=data(7,ii-1)-data(7,ii-2);
        dy=data(8,ii-1)-data(8,ii-2);
        cog=atan2d(dy,dx);
        D=[0 0 -1*tand(cog) 1 0 0];
        d=0;
        K=K-D'/(D*D')*(D*x_k_minus-d)/(r_k'/S*r_k)*(r_k'/S);
    end
    
    %Update state vector
    x_k_plus=x_k_minus+K*r_k;
    
    %if ii<=6
    %    field="xp"+ii;
    %    rts.(field)=x_k_plus;
    %else
    %    rts.xp1=rts.xp2;
    %    rts.xp2=rts.xp3;
    %    rts.xp3=rts.xp4;
    %    rts.xp4=rts.xp5;
    %    rts.xp5=rts.xp6;
    %    rts.xp6=x_k_plus;
    %end
    
    %Velocity check
    %if ii>1
    %    delta_x=x_k_plus(1,1)-data(7,ii-1);
    %    delta_y=x_k_plus(2,1)-data(8,ii-1);
    %    velocity=hypot(delta_x,delta_y)/dt;
    %    data(24,ii)=velocity;
    %    if delta_x>20
    %        x_k_plus(1,1)=x_k_plus(1,1)-delta_x+20;
    %    elseif delta_x<-20
    %        x_k_plus(1,1)=x_k_plus(1,1)-delta_x-20;
    %    end
    %    if delta_y>20
    %        x_k_plus(2,1)=x_k_plus(2,1)-delta_y+20;
    %    elseif delta_y<-20
    %        x_k_plus(2,1)=x_k_plus(2,1)-delta_y-20;
    %    end
    %end
    
    %Estimate projection with non-linear constraint (stw)
    %g_x=hypot((x_k_minus(3,1)-x_k_minus(5,1)),(x_k_minus(4,1)-x_k_minus(6,1)));
    %dg_dx=(x_k_minus(3,1)-x_k_minus(5,1))/g_x;
    %dg_dy=(x_k_minus(4,1)-x_k_minus(6,1))/g_x;
    %D=[0 0 dg_dx dg_dy -1*dg_dx -1*dg_dy];
    %d=stw-g_x+(D*x_k_minus);
    %x_k_plus=x_k_plus-D'/(D*D')*(D*x_k_plus-d);
    
    %Estimate projection with non-linear constraint (range in x)
    %g_x=x_k_minus(1,1)/cosd(z_k(2,1));
    %dg_dx=1/cosd(z_k(2,1));
    %D=[dg_dx 0 0 0 0 0];
    %d=z_k(1,1)-g_x+(D*x_k_minus);
    %x_k_plus=x_k_plus-D'/(D*D')*(D*x_k_plus-d);
    
    %Estimate projection with non-linear constraint (range in y)
    %g_x=x_k_minus(2,1)/sind(z_k(2,1));
    %dg_dy=1/sind(z_k(2,1));
    %D=[0 dg_dy 0 0 0 0];
    %d=z_k(1,1)-g_x+(D*x_k_minus);
    %x_k_plus=x_k_plus-D'/(D*D')*(D*x_k_plus-d);
    
    %Estimate projection with non-linear constraint(azi)
    %g_x=x_k_minus(2,1)/x_k_minus(1,1);
    %dg_dx=-1*x_k_minus(2,1)/(x_k_minus(1,1)^2);
    %dg_dy=1/x_k_minus(1,1);
    %D=[dg_dx dg_dy 0 0 0 0];
    %d=tand(z_k(2,1))-g_x+(D*x_k_minus);
    %x_k_plus=x_k_plus-D'/(D*D')*(D*x_k_plus-d);
    
    %Estimate projection with linear constraint (cog)
    %if ii>2
    %    dx=data(7,ii-1)-data(7,ii-2);
    %    dy=data(8,ii-1)-data(8,ii-2);
    %    cog=atan2d(dy,dx);
    %    D=[0 0 -1*tand(cog) 1 0 0];
    %    d=0;
    %    x_k_plus=x_k_plus-D'/(D*D')*(D*x_k_plus-d);
    %end
   
    %Simulate affect of GPS reset
    %Did not help stabilize the norm of P
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
    
    %if ii<=6
    %    field="pp"+ii;
    %    rts.(field)=P_plus;
    %else
    %    rts.pp1=rts.pp2;
    %    rts.pp2=rts.pp3;
    %    rts.pp3=rts.pp4;
    %    rts.pp4=rts.pp5;
    %    rts.pp5=rts.pp6;
    %    rts.pp6=P_plus;
    %end
    
    %Periodically reset velocities to the mean
    %Seemed to help stabilize error when one course was used.  Was not
    %useful with multiple courses
    %remain=mod(ii,300);
    %if remain==0
    %    x_k_plus(3,1)=mean(data(9,:));
    %    x_k_plus(4,1)=mean(data(10,:));
    %    x_k_plus(5,1)=mean(data(11,:));
    %    x_k_plus(6,1)=mean(data(12,:));
    %end
    
    %RTS smoother
    %if ii>=6
    %    x_nn=x_k_plus;
    %    P_nn=P_plus;
    %    for nn=5:-1:1
    %        pp_ind="pp"+nn;
    %        pm_ind="pm"+(nn+1);
    %        xm_ind="xm"+(nn+1);
    %        xp_ind="xp"+nn;
    %        
    %        K_nn=rts.(pp_ind)*F'/rts.(pm_ind);
    %        P_nn=rts.(pp_ind)-K_nn*(rts.(pm_ind)-P_nn)*K_nn';
    %        x_nn=rts.(xp_ind)+K_nn*(x_nn-rts.(xm_ind));
    %    end
    %    rts_out(:,ii-5)=x_nn;
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
    data(19,ii)=azi-atan2d(x_k_plus(2,1),x_k_plus(1,1));
    data(20,ii)=S(1,1);
    data(21,ii)=S(2,1);
    data(22,ii)=r_k(1,1);
    data(23,ii)=r_k(2,1);
end

avg(1,1)=mean(data(11,1:2000));
avg(2,1)=mean(data(12,1:2000));
avg(3,1)=atan2d(mean(data(12,1:2000)),mean(data(11,1:2000)));
avg(4,1)=mean(data(11,2001:steps));
avg(5,1)=mean(data(12,2001:steps));
avg(6,1)=atan2d(mean(data(12,2001:steps)),mean(data(11,2001:steps)));

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

%figure
%plot(rts_out(1,:),rts_out(2,:),data(15,6:steps),data(16,6:steps))
%title('Acutal Position vs RTS Output Position')
%xlabel('x position (m)')
%ylabel('y position (m)')
%legend('RTS Output','Actual Position','location','southeast')

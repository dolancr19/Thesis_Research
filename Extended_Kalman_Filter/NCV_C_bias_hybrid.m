function [data,epsilon,epsilon_v,mu, K_out, R_out, P_plus_out, P_minus_out]= NCV_C_bias_hybrid(freq, steps, x0_pos, y0_pos, stw0, hdg0, track, current,v_k_full,track0,stw_bias)
%% Prepare workspace (Uncomment this section to run one iteration for troubleshooting)
% clear variables;
% clc
% 
% % Define contants
% 
% freq=1; %Number of cycles per second
% steps=3600*freq; %# of filter steps desired
% track0=build_track(freq,steps); %build simulated values for speed through the water and heading
% stw_bias=.1;
% track.hdg=track0.hdg;
% track.stw=track0.stw+stw_bias;
% % 
% % v_k_full=randn(4,steps);
% set0=.25;
% drift0=90-60;
% variation=0; % zero for constant current, non-zero defines band for oscillating current
% current=build_current(steps, set0, drift0, variation);
% % rng(1);
% % s=rng;
% % % data_full=zeros(27,steps,iterations); %matrix for post-processing data
% % % epsilon_full=zeros(iterations,steps);
% % epsilon_bar=zeros(1,steps);
% % % epsilon_v_full=zeros(iterations,steps);
% % epsilon_v_bar=zeros(1,steps);
% % % mu_full=zeros(6,steps,iterations);
% % mu_bar=zeros(6,steps);
% % 
% %Define boundary conditions
% x0_pos=0;
% y0_pos=0;
% stw0=.50;
% hdg0=90-30;
%% Initilize constants, vectors and matrices
dt=1/freq; %Define time step size
data=zeros(29,steps); %matrix for post-processing data

%Initialize state vectors
x0_vel=stw0*cosd(hdg0);
y0_vel=stw0*sind(hdg0);

e0_current=0;%current.set(1)*cosd(current.drift(1));
n0_current=0;%current.set(1)*sind(current.drift(1));
%Initialize state vectors
x_k_minus=[x0_pos;x0_vel;y0_pos;y0_vel;e0_current;n0_current;1]; 
x_k_plus=x_k_minus;

%Define Q, the covariance matrix for noise associated with the state vector
var_Qp=.575; % position .575 last baseline
var_Qw=.025; % water referenced velocity, .035 is baseline
var_Qc=.01; % current velocity, .01 is baseline
var_Qb=.001; % .005 is baseline
Q=diag([var_Qp^2,var_Qw^2,var_Qp^2,var_Qw^2,var_Qc^2,var_Qc^2,var_Qb^2]);

%Compute R, the covariance matrix associated with measurement error
var_Rr_r=1^2; % range 1 last baseline
var_Ra_r=deg2rad(.55^2); % azimuth .55 last baseline
var_Rs_r=.5^2; % speed through the water, .5 is baseline
var_Ra1_r=deg2rad(.55^2); % heading .55 last baseline
%radians

R=zeros(4);
R(3,3)=var_Rs_r;
R(4,4)=var_Ra1_r;

%Initialize error covariance matrices
var_Pr=0.5^2; % position
var_Ps=.05^2; % velocity
var_Pb=.5^2; % bias
P_minus=diag([var_Pr,var_Ps,var_Pr,var_Ps,var_Ps,var_Ps, var_Pb]);
P_plus=P_minus;

%Initialize measurement vector
range_x=x0_pos;
range_y=y0_pos;

%Initialize state transition matrix
F=[1 dt 0 0 dt 0 0 ;0 1 0 0 0 0 0 ;0 0 1 dt 0 dt 0 ;0 0 0 1 0 0 0 ;0 0 0 0 1 0 0 ;0 0 0 0 0 1 0 ;0 0 0 0 0 0 1];

%Initialize vectors for consistency checks
epsilon=zeros(1,steps);
epsilon_v=zeros(1,steps);
mu=zeros(length(x_k_plus),steps);
K_out=zeros(length(x_k_plus),steps);
R_out=zeros(4,4,steps);
P_minus_out=zeros(length(x_k_plus),steps);
P_plus_out=zeros(length(x_k_plus),steps);

%% Execute filter
for ii=1:steps
    %Calculate state estimate and error covariance matrix for next step
    x_k_minus=F*x_k_plus;
    P_minus=F*P_plus*F.'+Q;
    P_minus_out(:,ii)=diag(P_minus);    
    data(1:7,ii)=x_k_minus;
              
    %Calculate noisy measurement vector
    stw0=track0.stw(ii);
    hdg0=track0.hdg(ii);
    stw0_x=stw0*cosd(hdg0);
    stw0_y=stw0*sind(hdg0);
    
    stw=track.stw(ii);
    hdg=track.hdg(ii);
    
    bias_ratio=stw/stw0;
    
    x_vel=stw0_x+current.set(ii)*cosd(current.drift(ii));
    y_vel=stw0_y+current.set(ii)*sind(current.drift(ii));
    range_x=range_x+x_vel*dt;
    range_y=range_y+y_vel*dt;
    range=hypot(range_x,range_y);
    azi=atan2d(range_y,range_x);
    azi_r=deg2rad(azi);
    
    x_act=[range_x;stw0_x;range_y;stw0_y;current.set(ii)*cosd(current.drift(ii));current.set(ii)*sind(current.drift(ii));bias_ratio];
 
    %Executes the filter at the specified frequency with one second
    %measurement updates.  Simulates a 1 Hz iUSBL transmission.
    if ii==1 || mod(ii,1)==0
        % Calculate converted measurement bias
        mu_t=[range*cos(azi_r)*((exp(-.5*var_Ra_r))-1);range*sin(azi_r)*((exp(-.5*var_Ra_r))-1);0;0];

        %Calculate converted measurement covariances
        R(1,1)=(range^2*exp(-1*var_Ra_r)*((cos(azi_r)^2)*(cosh(var_Ra_r)-1)+(sin(azi_r)^2*sinh(var_Ra_r))))+(var_Rr_r*exp(-1*var_Ra_r)*((cos(azi_r)^2)*(cosh(var_Ra_r))+(sin(azi_r)^2*sinh(var_Ra_r))));
        R(2,2)=(range^2*exp(-1*var_Ra_r)*((sin(azi_r)^2)*(cosh(var_Ra_r)-1)+(cos(azi_r)^2*sinh(var_Ra_r))))+(var_Rr_r*exp(-1*var_Ra_r)*((sin(azi_r)^2)*(cosh(var_Ra_r))+(cos(azi_r)^2*sinh(var_Ra_r))));   

        R(1,2)=sin(azi_r)*cos(azi_r)*exp(-2*var_Ra_r)*(var_Rr_r+range^2*(1-exp(var_Ra_r)));
        R(2,1)=R(1,2);

        v_k=randn(4,1);
        [V_R,D_R]=eig(R);
        [~,ind] = sort(diag(D_R),'descend');
        D_Rs = D_R(ind,ind);
        V_Rs = V_R(:,ind);
        d=diag(D_Rs);
        d=sqrt(d);
        v_k=V_Rs*(d.*v_k); %varible noise for Monte Carlo testing
    %     v_k=V_R*(d.*v_k_full(:,ii)); %same noise for bias testing
        z_k=[range_x;range_y;stw;hdg]+v_k-mu_t;

        %Calculate measurement vector, both linear and non-linear elements
        h=[x_k_minus(1,1);x_k_minus(3,1);x_k_minus(7,1)*hypot(x_k_minus(2,1),x_k_minus(4,1));atan2d(x_k_minus(4,1),x_k_minus(2,1))];

        %Calculate measurement mapping matrix
        stw_denom=x_k_minus(7,1)/hypot(x_k_minus(2,1),x_k_minus(4,1));
        stw_dx=x_k_minus(2,1)*stw_denom;
        stw_dy=x_k_minus(4,1)*stw_denom;

        hdg_dx=-1*x_k_minus(4,1)/(x_k_minus(2,1)^2+x_k_minus(4,1)^2);
        hdg_dy=1/(x_k_minus(2,1)+(x_k_minus(4,1)^2/x_k_minus(2,1)));

        H=[1 0 0 0 0 0 0;0 0 1 0 0 0 0;0 stw_dx 0 stw_dy 0 0 1;0 hdg_dx 0 hdg_dy 0 0 0];
        
        %Calculate innovation
        r_k(1,1)=z_k(1,1)-h(1,1);
        r_k(2,1)=z_k(2,1)-h(2,1);
        r_k(3,1)=z_k(3,1)-h(3,1);
        r_k(4,1)=angdiff(deg2rad(h(4,1)),deg2rad(z_k(4,1)));

        %Calculate innovation covariance
        S=H*P_minus*H.'+R;

        %Update Kalman Gain after estimate vector calculation
        K=(P_minus*H.')/S;

        %Update state vector
        x_k_plus=x_k_minus+K*r_k;

        %Update error covariance matrix with new Kalman Gain
        %Joseph stabilized equation
        P_plus=(eye(length(x_k_plus))-K*H)*P_minus*((eye(length(x_k_plus))-K*H).')+K*R*K.';
    else
       x_k_plus=x_k_minus;
       P_plus=P_minus;
    end  
    
    
    P_plus_out(:,ii)=diag(P_plus);
    
    %Calculate the normalized estimation error squared (NEES)
    x_err=x_act-x_k_plus;
    epsilon(ii)=x_err.'*(P_plus\x_err);
    
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

end
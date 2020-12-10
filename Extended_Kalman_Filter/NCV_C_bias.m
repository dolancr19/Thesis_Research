function [data,epsilon,epsilon_v,mu, K_out, R_out, P_plus_out, P_minus_out]= NCV_C_bias(freq, steps, x0_pos, y0_pos, stw0, hdg0, track, current,v_k_full,track0,stw_bias)
%% Prepare workspace
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
%%
dt=1/freq; %Define time step size
data=zeros(29,steps); %matrix for post-processing data

% bias_mag=0;
% last=1;
x0_vel=stw0*cosd(hdg0);
y0_vel=stw0*sind(hdg0);

e0_current=0;%current.set(1)*cosd(current.drift(1));
n0_current=0;%current.set(1)*sind(current.drift(1));
%Initialize state vectors
x_k_minus=[x0_pos;x0_vel;y0_pos;y0_vel;e0_current;n0_current;1]; 
x_k_plus=x_k_minus;
% x_k_plus_RA=x_k_plus(1:4,1);

%Define Q, the covariance matrix for noise associated with the state vector
var_Qp=.575; % position .575 last baseline
var_Qw=.025; % water referenced velocity, .035 is baseline
var_Qc=.01; % current velocity, .01 is baseline
var_Qb=.01; % 1.5 * bias
Q=diag([var_Qp^2,var_Qw^2,var_Qp^2,var_Qw^2,var_Qc^2,var_Qc^2,var_Qb^2]);
% Q_RA=diag([var_Qp^2,var_Qw^2,var_Qp^2,var_Qw^2]);%,var_Qc^2,var_Qc^2]);

%Compute Q with white noise acceleration plus process noise for the current
% var_Q=.01^2;
% var_Qc=.1^2;
% Qc=diag([var_Qc var_Qc]); 
% %var_Q=.0001^2;
% G=[.5*(dt^2); dt;.5*(dt^2);dt];
% Q=G*var_Q*G';
% Q=[Q zeros(4,2);zeros(2,4) Qc];

%Compute R, the covariance matrix associated with measurement error
var_Rr_r=1^2; % range 1 last baseline
var_Ra_r=deg2rad(.55^2); % azimuth .55 last baseline
var_Rs_r=.5^2; % speed through the water, .5 is baseline
var_Ra1_r=deg2rad(.55^2); % heading .55 last baseline
%radians

R=zeros(4); %normal
% R=zeros(2); %remove simulated iUSBL measurements for bias testing

%Initialize error covariance matrices
var_Pr=0.5^2; % position
var_Ps=.05^2; % velocity
var_Pb=.01^2; % bias
P_minus=diag([var_Pr,var_Ps,var_Pr,var_Ps,var_Ps,var_Ps, var_Pb]);
P_plus=P_minus;

% P_minus_RA=diag([var_Pr,var_Ps,var_Pr,var_Ps]);%,var_Ps,var_Ps]);
% P_plus_RA=P_minus_RA;

% b=[0;0;0;0];
% b_cov=zeros(4);
%Initialize measurement vector
range_x=x0_pos;
range_y=y0_pos;

% Using water referenced speed in state vector
% F=[1 dt 0 0 dt 0 0 0;0 1 0 0 0 0 0 0;0 0 1 dt 0 dt 0 0;0 0 0 1 0 0 0 0;0 0 0 0 1 0 0 0;0 0 0 0 0 1 0 0;0 0 0 0 0 0 1 0;0 0 0 0 0 0 0 1]; % state transition matrix
% H=[1 0 0 0 0 0 0 0;0 1 0 0 0 0 1 0;0 0 1 0 0 0 0 0;0 0 0 1 0 0 0 1]; 

%using ground referenced speed in state vector
F=[1 dt 0 0 dt 0 0 ;0 1 0 0 0 0 0 ;0 0 1 dt 0 dt 0 ;0 0 0 1 0 0 0 ;0 0 0 0 1 0 0 ;0 0 0 0 0 1 0 ;0 0 0 0 0 0 1]; % state transition matrix
H=[1 0 0 0 0 0 0;0 1 0 0 0 0 0;0 0 1 0 0 0 0;0 0 0 1 0 0 0]; 

%Initialize vectors for consistency checks
epsilon=zeros(1,steps);
epsilon_v=zeros(1,steps);
mu=zeros(length(x_k_plus),steps);
K_out=zeros(length(x_k_plus),steps);
R_out=zeros(4,4,steps);
P_minus_out=zeros(length(x_k_plus),steps);
P_plus_out=zeros(length(x_k_plus),steps);
% b_out=zeros(4,steps);
%% Execute filter
for ii=1:steps
    %Calculate state estimate and error covariance matrix for next step
    x_k_minus=F*x_k_plus;
    bias_noise=randn(1,1)*.0075;%var_Qb;
    x_k_minus(7,1)=x_k_minus(7,1);%+bias_noise;
    P_minus=F*P_plus*F.'+Q;
    P_minus_out(:,ii)=diag(P_minus);    
    data(1:7,ii)=x_k_minus;
              
    %Calculate noisy measurement vector
    stw0=track0.stw(ii);
    hdg0=track0.hdg(ii);
    stw0_x=stw0*cosd(hdg0);
    stw0_y=stw0*sind(hdg0);
    
    
%     if mod(ii,60)==0
%         
%        bias_mag=mean(hypot(data(7,last:ii), data(8,last:ii)));
% %      bias_mag=hypot(x_k_minus(7,1), x_k_minus(8,1));   
%      var_Qb = 9.5*bias_mag^2 + -0.75*bias_mag + 0.12;
%         Q(8,8)=var_Qb^2;
%         last=ii;
%        
%    end
%     bias_mag=hypot(x_k_minus(7,1), x_k_minus(8,1));
%  bias_mag_out(ii)=bias_mag;
    stw=track.stw(ii);%-bias_mag;
    hdg=track.hdg(ii);
    
    bias_ratio=stw0/stw;
    
    var_Rs_rm=((var_Qb^2)*var_Rs_r)+((var_Qb^2)*stw^2)+(var_Rs_r*(x_k_minus(7,1)^2));
    stw=stw*x_k_minus(7,1);
    
    stw_x=stw*cosd(hdg);
    stw_y=stw*sind(hdg);
    
 
    x_vel=stw0_x+current.set(ii)*cosd(current.drift(ii));
    y_vel=stw0_y+current.set(ii)*sind(current.drift(ii));
    range_x=range_x+x_vel*dt;
    range_y=range_y+y_vel*dt;
    range=hypot(range_x,range_y);
    azi=atan2d(range_y,range_x);
    
    azi_r=deg2rad(azi);
    hdg_r=deg2rad(hdg);
    
    
        
    
    x_act=[range_x;stw0_x;range_y;stw0_y;current.set(ii)*cosd(current.drift(ii));current.set(ii)*sind(current.drift(ii));bias_ratio];
 
        %Executes the filter at the specified frequency with one second
    %measurement updates.  Simulates a 1 Hz iUSBL transmission.
    if ii==1 || mod(ii,1)==0
    % Calculate converted measurement bias
    %Normal
    mu_t=[range*cos(azi_r)*((exp(-.5*var_Ra_r))-1);stw*cos(hdg_r)*((exp(-.5*var_Ra1_r))-1);range*sin(azi_r)*((exp(-.5*var_Ra_r))-1);stw*sin(hdg_r)*((exp(-.5*var_Ra1_r))-1)];
    
    R(1,1)=(range^2*exp(-1*var_Ra_r)*((cos(azi_r)^2)*(cosh(var_Ra_r)-1)+(sin(azi_r)^2*sinh(var_Ra_r))))+(var_Rr_r*exp(-1*var_Ra_r)*((cos(azi_r)^2)*(cosh(var_Ra_r))+(sin(azi_r)^2*sinh(var_Ra_r))));
    R(2,2)=(stw^2*exp(-1*var_Ra1_r)*((cos(hdg_r)^2)*(cosh(var_Ra1_r)-1)+(sin(hdg_r)^2*sinh(var_Ra1_r))))+(var_Rs_rm*exp(-1*var_Ra1_r)*((cos(hdg_r)^2)*(cosh(var_Ra1_r))+(sin(hdg_r)^2*sinh(var_Ra1_r))));
    R(3,3)=(range^2*exp(-1*var_Ra_r)*((sin(azi_r)^2)*(cosh(var_Ra_r)-1)+(cos(azi_r)^2*sinh(var_Ra_r))))+(var_Rr_r*exp(-1*var_Ra_r)*((sin(azi_r)^2)*(cosh(var_Ra_r))+(cos(azi_r)^2*sinh(var_Ra_r))));
    R(4,4)=(stw^2*exp(-1*var_Ra1_r)*((sin(hdg_r)^2)*(cosh(var_Ra1_r)-1)+(cos(hdg_r)^2*sinh(var_Ra1_r))))+(var_Rs_rm*exp(-1*var_Ra1_r)*((sin(hdg_r)^2)*(cosh(var_Ra1_r))+(cos(hdg_r)^2*sinh(var_Ra1_r))));

    R(1,3)=sin(azi_r)*cos(azi_r)*exp(-2*var_Ra_r)*(var_Rr_r+range^2*(1-exp(var_Ra_r)));
    R(3,1)=R(1,3);
    
    R(2,4)=sin(hdg_r)*cos(hdg_r)*exp(-2*var_Ra1_r)*(var_Rs_rm+stw^2*(1-exp(var_Ra1_r)));
    R(4,2)=R(2,4);

       
%     R=R+b_cov;

%     %Remove simulated iUSBL measurements for bias testing
%     mu_t=[stw*cos(hdg_r)*((exp(-.5*var_Ra1_r))-1);stw*sin(hdg_r)*((exp(-.5*var_Ra1_r))-1)];
%     
%     R(1,1)=(stw^2*exp(-1*var_Ra1_r)*((cos(hdg_r)^2)*(cosh(var_Ra1_r)-1)+(sin(hdg_r)^2*sinh(var_Ra1_r))))+(var_Rs_r*exp(-1*var_Ra1_r)*((cos(hdg_r)^2)*(cosh(var_Ra1_r))+(sin(hdg_r)^2*sinh(var_Ra1_r))));
%     R(2,2)=(stw^2*exp(-1*var_Ra1_r)*((sin(hdg_r)^2)*(cosh(var_Ra1_r)-1)+(cos(hdg_r)^2*sinh(var_Ra1_r))))+(var_Rs_r*exp(-1*var_Ra1_r)*((sin(hdg_r)^2)*(cosh(var_Ra1_r))+(cos(hdg_r)^2*sinh(var_Ra1_r))));
% 
%     R(1,2)=sin(hdg_r)*cos(hdg_r)*exp(-2*var_Ra1_r)*(var_Rs_r+stw^2*(1-exp(var_Ra1_r)));
%     R(2,1)=R(1,2);

%     % Linearized covariance matrix
%     mu_t=0;
%     
%     R(1,1)=(range^2*var_Ra_r*(sin(azi_r)^2))+(var_Rr_r*(cos(azi_r)^2));
%     R(2,2)=(stw^2*var_Ra1_r*(sin(hdg_r)^2))+(var_Rs_r*(cos(hdg_r)^2));
%     R(3,3)=(range^2*var_Ra_r*(cos(azi_r)^2))+(var_Rr_r*(sin(azi_r)^2));
%     R(4,4)=(stw^2*var_Ra1_r*(cos(hdg_r)^2))+(var_Rs_r*(sin(hdg_r)^2));
%     
%     R(1,3)=(var_Rr_r-range^2*var_Ra_r)*sin(azi_r)*cos(azi_r);
%     R(3,1)=R(1,3);
%     
%     R(2,4)=(var_Rs_r-stw^2*var_Ra1_r)*sin(hdg_r)*cos(hdg_r);
%     R(4,2)=R(2,4);
    
     v_k=randn(4,1); %normal
%     v_k=randn(2,1); %remove simulated iUSBL measurements for bias testing
    [V_R,D_R]=eig(R);
%         [~,ind] = sort(diag(D_R),'descend');
%         D_Rs = D_R(ind,ind);
%         V_Rs = V_R(:,ind);
    d=diag(D_R);
    %d=real(sqrt(d));
    d=sqrt(d);
    v_k=V_R*(d.*v_k); %varible noise for Monte Carlo testing
%     v_k=V_R*(d.*v_k_full(:,ii)); %same noise for bias testing
    z_k=[range_x;stw_x;range_y;stw_y]+v_k-mu_t;%+b; %normal
%     z_k=[stw_x;stw_y]+v_k-mu_t; %remove simulated iUSBL measurements for bias testing

    %Calculate innovation
    r_k=z_k-H*x_k_minus;

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
    
%     [x_k_plus_RA,P_plus_RA]= NCV_C_RA(freq,x_k_plus_RA,P_plus_RA, Q_RA,v_k,range_x,range_y);
%     
%     sog_all=[x_k_plus(2,1);x_k_plus(4,1)]+[x_k_plus(5,1);x_k_plus(6,1)];
%     sog_RA=[x_k_plus_RA(2,1);x_k_plus_RA(4,1)];
%     
%     sog_diff=sog_RA-sog_all;
%     
%     b(2,1)=sog_diff(1,1);
%     b(4,1)=sog_diff(2,1);
%     b_out(:,ii)=b;
%     
%     sog_cov=[P_plus(2,2) P_plus(2,4);P_plus(4,2) P_plus(4,4)]+[P_plus(5,5) P_plus(5,6);P_plus(6,5) P_plus(6,6)]+[P_plus_RA(2,2) P_plus_RA(2,4);P_plus_RA(4,2) P_plus_RA(4,4)];
%     
%     b_cov(2,2)=sog_cov(1,1);
%     b_cov(4,4)=sog_cov(2,2);
%     b_cov(2,4)=sog_cov(1,2);
%     b_cov(4,2)=sog_cov(2,1);
    
    %P_plus=P_plus+b_cov;
    
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
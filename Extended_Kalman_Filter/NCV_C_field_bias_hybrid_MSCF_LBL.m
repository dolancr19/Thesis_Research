%function [data,epsilon_v]= NCV_field(NAVDR, SOURCEXY, ACOUSTICRB, e_gps_interp, n_gps_interp) %MOOS data
function [data,epsilon_v,rho_bar,P_minus_out,z_k_out,F_out,K_out2,K_out4,acoustic_start,acoustic_end]=NCV_C_field_bias_hybrid_MSCF_LBL(NAVDR,LBL_NAV,SOURCE_XY,NAVXY,first,last) %PF data
    steps=last;%length(NAVDR.TIME);
    data=zeros(21,steps); %matrix for post-processing data

    %Define boundary conditions

    acoustic_count=1;

    while NAVDR.TIME(first)-LBL_NAV.time(acoustic_count)>0
        acoustic_count=acoustic_count+1;
    end
    acoustic_count=acoustic_count-1;
    acoustic_start=acoustic_count;
    
    e0_pos=LBL_NAV.xs(acoustic_count);
    n0_pos=LBL_NAV.ys(acoustic_count);

    stw0=NAVDR.NAV_SPEED(first);
    hdg0=90-NAVDR.NAV_HEADING(first);

    e0_vel=stw0*cosd(hdg0);
    n0_vel=stw0*sind(hdg0);


% acoustic_bearing=90-(NAVDR.NAV_HEADING(1)-rad2deg(NAVDR.Bearing(1)));
% if acoustic_bearing>180
%     acoustic_bearing=acoustic_bearing-360;
% elseif acoustic_bearing<-180
%     acoustic_bearing=acoustic_bearing+360;
% end
% source_x=Interp.Range(1)*cosd(acoustic_bearing);
% source_y=Interp.Range(1)*sind(acoustic_bearing);
% 
% e0_pos=Interp.e_gps(1)-source_x;
% n0_pos=Interp.n_gps(1)-source_y;
% 
% stw0=Interp.NAV_SPEED(1);
% hdg0=90-Interp.NAV_HEADING(1);
% 
% e0_vel=stw0*cosd(hdg0);
% n0_vel=stw0*sind(hdg0);

    %Initialize state and measurement vectors
    x_k_minus=[e0_pos;e0_vel;n0_pos;n0_vel;0;0;1.63];
    data(1:7,first)=x_k_minus;
    x_k_plus=x_k_minus;
    data(8:14,first)=x_k_plus;
    z_k_out=zeros(4,steps);
    z_k_out(:,first)=[e0_pos,n0_pos,stw0,hdg0];
% r_k=zeros(4,1);

    %Define Q, the covariance matrix for noise associated with the state vector
    var_Qp=1; % position 5, 1 quokka
    var_Qw=.5; % water referenced velocity .1
    var_Qc=.05; % current velocity .01, .1 quokka
    var_Qb=.001; % bias .001
    Q=diag([var_Qp^2,var_Qw^2,var_Qp^2,var_Qw^2,var_Qc^2,var_Qc^2,var_Qb^2]);

    %Define R, the covariance matrix associated with measurement error, value will be
    %calculated during each step of the filter
    var_Rr_r=10^2; % range 15, 10 quokka
%     var_Ra_r=deg2rad(4^2); % azimuth 2
    var_Rs_r=.5^2; % speed through the water .5
    var_Ra1_r=deg2rad(2^2); % heading 1, 3 quokka



    %Initialize error covariance matrices
    var_Pr=10^2; % position
    var_Ps=10^2; % velocity
    var_Pb=.5^2; % bias
    P_minus=diag([var_Pr,var_Ps,var_Pr,var_Ps,var_Ps,var_Ps,var_Pb]);
    P_plus=P_minus;

    sigma_k=2;% .9 quokka

   

    %Initialize vectors for consistency checks
    %epsilon=zeros(1,steps);
    epsilon_v=zeros(1,steps);
    %mu=zeros(length(x_k_plus),steps);
    K_out2=zeros(length(x_k_minus),2,steps);
    K_out4=zeros(length(x_k_minus),4,steps);
    % S_out=zeros(length(R),length(R),steps);
    P_minus_out=zeros(length(x_k_minus),length(x_k_minus),steps);
    F_out=zeros(length(x_k_minus),length(x_k_minus),steps);
%% Execute filter
    for ii=first+1:steps
        dt=NAVDR.TIME(ii)-NAVDR.TIME(ii-1);
        F=[1 dt 0 0 dt 0 0 ;0 1 0 0 0 0 0 ;0 0 1 dt 0 dt 0 ;0 0 0 1 0 0 0 ;0 0 0 0 1 0 0 ;0 0 0 0 0 1 0 ;0 0 0 0 0 0 1];
%         F_out(:,:,ii)=F;

        
        
        %Calculate state estimate and error covariance matrix for next step
        x_k_minus=F*x_k_plus;
        P_minus=F*P_plus*F.'+Q;
%         P_minus_out(:,:,ii)=P_minus;    
        data(1:7,ii)=x_k_minus;
        if NAVDR.NAV_SPEED(ii)==0
            x_k_plus=x_k_minus;
            P_plus=P_minus;
        else
        
            
            if NAVDR.TIME(ii)-LBL_NAV.time(acoustic_count)<0
                [x_k_plus,P_plus]=no_range();
            else
                if isnan(LBL_NAV.xs(acoustic_count)) 
                    %Calculate measurement vec and vector
                    [x_k_plus,P_plus]=no_range();
                else
                    [x_k_plus,P_plus]=full_measurement();
                end
                acoustic_count=acoustic_count+1;
                if acoustic_count>length(LBL_NAV.time)
                    acoustic_count=acoustic_count-1;
                end
            end
        end
        data(8:14,ii)=x_k_plus;
        data(15,ii)=hypot(x_k_plus(2,1),x_k_plus(4,1));
        data(16,ii)=atan2d(x_k_plus(4,1),x_k_plus(2,1));
        data(17,ii)=norm(P_plus);
    end

    rho_data=data(18:21,:);

    trim_r=zeros(1,length(rho_data(1,:)));
    for rr=1:length(trim_r)
        if rho_data(1,rr)==0
            trim_r(rr)=0;
        else
            trim_r(rr)=1;
        end
    end
%     trim_r=logical(trim_r);
%     rho_data_trim=rho_data(:,trim_r);
% 
%     sum_kj=zeros(4,1);
%     sum_kk=zeros(4,1);
%     sum_jj=zeros(4,1);
%     for jj=1:length(rho_data_trim)-1
%         sum_kj(1)=sum_kj(1)+(rho_data_trim(1,jj)*rho_data_trim(1,jj+1));
%         sum_kj(2)=sum_kj(2)+(rho_data_trim(2,jj)*rho_data_trim(2,jj+1));
%         sum_kj(3)=sum_kj(3)+(rho_data_trim(3,jj)*rho_data_trim(3,jj+1));
%         sum_kj(4)=sum_kj(4)+(rho_data_trim(4,jj)*rho_data_trim(4,jj+1));
%         sum_kk(1)=sum_kk(1)+(rho_data_trim(1,jj))^2;
%         sum_kk(2)=sum_kk(2)+(rho_data_trim(2,jj))^2;
%         sum_kk(3)=sum_kk(3)+(rho_data_trim(3,jj))^2;
%         sum_kk(4)=sum_kk(4)+(rho_data_trim(4,jj))^2;
%         sum_jj(1)=sum_jj(1)+(rho_data_trim(1,jj+1))^2;
%         sum_jj(2)=sum_jj(2)+(rho_data_trim(2,jj+1))^2;
%         sum_jj(3)=sum_jj(3)+(rho_data_trim(3,jj+1))^2;
%         sum_jj(4)=sum_jj(4)+(rho_data_trim(4,jj+1))^2;
%     end
    rho_bar=zeros(4,1);
%     rho_bar(1,1)=sum_kj(1)/sqrt(sum_kk(1)*sum_jj(1));
%     rho_bar(2,1)=sum_kj(2)/sqrt(sum_kk(2)*sum_jj(2));
%     rho_bar(3,1)=sum_kj(3)/sqrt(sum_kk(3)*sum_jj(3));
%     rho_bar(4,1)=sum_kj(4)/sqrt(sum_kk(4)*sum_jj(4));

    acoustic_end=acoustic_count;

    function [x_k_plus,P_plus]=no_range()
        stw=NAVDR.NAV_SPEED(ii);
        hdg=90-NAVDR.NAV_HEADING(ii);

        R=zeros(2);
        R(1,1)=var_Rs_r;
        R(2,2)=var_Ra1_r;

        R_in=inv(R);

        z_k=[stw;hdg];
        z_k_out(3:4,ii)=z_k;

        %Calculate measurement vector, both linear and non-linear elements
        h=[x_k_minus(7,1)*hypot(x_k_minus(2,1),x_k_minus(4,1));atan2d(x_k_minus(4,1),x_k_minus(2,1))];

        if hypot(x_k_minus(2,1),x_k_minus(4,1))==0
            stw_dx=0;
            stw_dy=0;
            hdg_dx=0;
            hdg_dy=0;
        else
            %Calculate measurement mapping matrix
            stw_denom=x_k_minus(7,1)/hypot(x_k_minus(2,1),x_k_minus(4,1));
            stw_dx=x_k_minus(2,1)*stw_denom;
            stw_dy=x_k_minus(4,1)*stw_denom;

            hdg_dx=-1*x_k_minus(4,1)/(x_k_minus(2,1)^2+x_k_minus(4,1)^2);
            hdg_dy=1/(x_k_minus(2,1)+(x_k_minus(4,1)^2/x_k_minus(2,1)));
        end

        H=[0 stw_dx 0 stw_dy 0 0 1;0 hdg_dx 0 hdg_dy 0 0 0];

        %Calculate innovation
        r_k=zeros(2,1);
        r_k(1,1)=z_k(1,1)-h(1,1);
        r_k(2,1)=angdiff(deg2rad(h(2,1)),deg2rad(z_k(2,1)));

        C_k=diag([exp(-1*R_in(1,1)*(r_k(1,1)^2)/(2*(sigma_k^2))),exp(-1*R_in(2,2)*(r_k(2,1)^2)/(2*(sigma_k^2)))]);

        K=inv((inv(P_minus)+H.'*C_k*R_in*H))*(H.'*C_k*R_in);
        if hypot(x_k_minus(2,1),x_k_minus(4,1))==0
            K=[0 0;.5 0;0 0;.5 0;0 0;0 0;0 0];
        end
        %Calculate innovation covariance
        S=H*P_minus*H.'+R;
    %     S_out(:,:,ii)=S;

        %Update Kalman Gain after estimate vector calculation
    %     K=(P_minus*H.')/S;
%         K_out2(:,:,ii)=K;

        %Update state vector
        x_k_plus=x_k_minus+K*r_k;

        %Update error covariance matrix with new Kalman Gain
        %Joseph stabilized equation
        P_plus=(eye(length(x_k_minus))-K*H)*P_minus*(eye(length(x_k_minus))-K*H)'+K*R*K';

        %Calculate the normalized innovation squared (NIS)
        epsilon_v(ii)=r_k.'*inv(S)*r_k;

        data(18:19,ii)=r_k;
    end

    function [x_k_plus,P_plus]=full_measurement()
        %Calculate measurement vec and vector
        
        stw=NAVDR.NAV_SPEED(ii);
        hdg=90-NAVDR.NAV_HEADING(ii);


        R=diag([var_Rr_r,var_Rr_r,var_Rs_r,var_Ra1_r]);

        R_inverse=inv(R);

        z_k=[LBL_NAV.xs(acoustic_count);LBL_NAV.ys(acoustic_count);stw;hdg];
        z_k_out(:,ii)=z_k;

        %Calculate measurement vector, both linear and non-linear elements
        h=[x_k_minus(1,1);x_k_minus(3,1);x_k_minus(7,1)*hypot(x_k_minus(2,1),x_k_minus(4,1));atan2d(x_k_minus(4,1),x_k_minus(2,1))];

    %     range_check(ii)=abs(hypot(h(1,1),h(2,1))-hypot(z_k(1,1),z_k(2,1)));
    %     if abs(hypot(h(1,1),h(2,1))-hypot(z_k(1,1),z_k(2,1)))>25
    %         [x_k_plus,P_plus]=no_range();
    %     else    

        %Calculate measurement mapping matrix
        stw_denom=x_k_minus(7,1)/hypot(x_k_minus(2,1),x_k_minus(4,1));
        stw_dx=x_k_minus(2,1)*stw_denom;
        stw_dy=x_k_minus(4,1)*stw_denom;

        hdg_dx=-1*x_k_minus(4,1)/(x_k_minus(2,1)^2+x_k_minus(4,1)^2);
        hdg_dy=1/(x_k_minus(2,1)+(x_k_minus(4,1)^2/x_k_minus(2,1)));

        H=[1 0 0 0 0 0 0;0 0 1 0 0 0 0;0 stw_dx 0 stw_dy 0 0 1;0 hdg_dx 0 hdg_dy 0 0 0];

        %Calculate innovation
    %     r_k=z_k-H*x_k_minus;
        r_k=zeros(4,1);
        r_k(1,1)=z_k(1,1)-h(1,1);
        r_k(2,1)=z_k(2,1)-h(2,1);
        r_k(3,1)=z_k(3,1)-h(3,1);
        r_k(4,1)=angdiff(deg2rad(h(4,1)),deg2rad(z_k(4,1)));

        C_k=diag([exp(-1*R_inverse(1,1)*(r_k(1,1)^2)/(2*(sigma_k^2))),exp(-1*R_inverse(2,2)*(r_k(2,1)^2)/(2*(sigma_k^2))),exp(-1*R_inverse(3,3)*(r_k(3,1)^2)/(2*(sigma_k^2))),exp(-1*R_inverse(4,4)*(r_k(4,1)^2)/(2*(sigma_k^2)))]);

        K=inv((inv(P_minus)+H.'*C_k*R_inverse*H))*(H.'*C_k*R_inverse);

        %Calculate innovation covariance
        S=H*P_minus*H.'+R;
    %     S_out(:,:,ii)=S;

        %Update Kalman Gain after estimate vector calculation
    %     K=(P_minus*H.')/S;
%         K_out4(:,:,ii)=K;

        %Update state vector
        x_k_plus=x_k_minus+K*r_k;

        %Update error covariance matrix with new Kalman Gain
        %Joseph stabilized equation
        P_plus=(eye(length(x_k_minus))-K*H)*P_minus*(eye(length(x_k_minus))-K*H)'+K*R*K';

        %Calculate the normalized innovation squared (NIS)
        epsilon_v(ii)=r_k.'*inv(S)*r_k;

        data(18:21,ii)=r_k;
    end
% end
end
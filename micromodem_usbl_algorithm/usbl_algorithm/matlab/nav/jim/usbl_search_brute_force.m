function [az_hat_deg, el_hat_deg, ximax, yimax] = usbl_search_brute_force(X, xyz_m, k0_1m, axvec, ayvec)

%axvec = (-1:0.01:+1)';
%ayvec = (-1:0.01:+1)';

Pval = zeros(length(axvec), length(ayvec));
ximax=NaN;
yimax=NaN;
maxPval = 0;
for xi=(1:length(axvec))
	ax = axvec(xi);
    for yi=(1:length(ayvec))
    	ay = ayvec(yi);
        az = sqrt(1-(ax^2)-(ay^2));
        k_vec_xi_yi = k0_1m * [ax; ay; az];
        S = transpose(exp(-sqrt(-1)*k_vec_xi_yi'*xyz_m));
        S_H = conj(transpose(S));
        Pval(xi,yi) = abs( S_H * X );
        if (Pval(xi,yi) > maxPval)
        	maxPval = Pval(xi,yi);
            ximax = xi;
            yimax = yi;
            %Smax = S;
            %k_vec_ximax_yimax = k_vec_xi_yi;
        end
    end
end

axhat = axvec(ximax);
ayhat = ayvec(yimax);
azhat2 = 1-(axhat^2)-(ayhat^2);
if (azhat2<0)
    fprintf('azhat2=%f; negative; zeroing\n', azhat2);
    azhat2=0;
end;
azhat = sqrt(azhat2);

az_hat_deg = (180/pi)*atan2(axhat,azhat);
el_hat_deg = (180/pi)*atan2(ayhat,azhat);
return;

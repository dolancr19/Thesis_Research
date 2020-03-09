function [r_out] = get_r(thetas, z_glider, r_glider)
load output
rss_compare = zeros(length(thetas),1);
r_compare = zeros(length(thetas),4);
for ii = 1:length(thetas)
    jj = 2;
    while rray_out(ii,jj) ~=0 
        if abs(rray_out(ii,jj)) < r_glider
            jj = jj+1;
        else
            jj = jj+1;
        end
    end
    jj=jj-1;
    r_compare1 = abs(r_glider-rray_out(ii,jj))/r_glider;
    r_compare2 = abs(r_glider-rray_out(ii,jj-1))/r_glider;
    if r_compare1 < r_compare2
        z_compare = abs(z_glider-zray_out(ii,jj))/z_glider;
        r_compare(ii,1) = rray_out(ii,jj);
        r_compare(ii,2) = zray_out(ii,jj);
        r_compare(ii,3) = tau_out(ii,jj);
        r_compare(ii,4) = thetavec_out(ii,jj);
        rss_compare(ii,1) = sqrt(r_compare1^2 + z_compare^2);
    else 
        z_compare = abs(z_glider-zray_out(ii,jj-1))/z_glider;
        r_compare(ii,1) = rray_out(ii,jj-1);
        r_compare(ii,2) = zray_out(ii,jj-1);
        r_compare(ii,3) = tau_out(ii,jj-1);
        r_compare(ii,4) = thetavec_out(ii,jj-1);
        rss_compare(ii,1) = sqrt(r_compare2^2 + z_compare^2);
    end
end
[~,ind] = sort(rss_compare, 'ascend');
r_compare_sort = r_compare(ind,:);
r_out = r_compare_sort(1,:);
end
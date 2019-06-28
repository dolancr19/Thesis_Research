clear variables
load output
thetas = 31;
z_glider = 2000;
tau_glider = 1.3716;
rss_compare = zeros(thetas,1);
r_compare = zeros(thetas,1);

for ii = 1:thetas
    jj = 2;
    while tau_out(ii,jj) >0 
        if tau_out(ii,jj) < tau_glider
            jj = jj+1;
        else
            jj = jj+1;
        end
    end
    jj=jj-1;
    
    tau_compare1 = abs(tau_glider-tau_out(ii,jj))/tau_glider;
    tau_compare2 = abs(tau_glider-tau_out(ii,jj-1))/tau_glider;
        
    if tau_compare1 < tau_compare2
        z_compare = abs(z_glider-zray_out(ii,jj))/z_glider;
        r_compare(ii,1) = rray_out(ii,jj);
        rss_compare(ii,1) = sqrt(tau_compare1^2 + z_compare^2);
    else 
        z_compare = abs(z_glider-zray_out(ii,jj-1))/z_glider;
        r_compare(ii,1) = rray_out(ii,jj-1);
        rss_compare(ii,1) = sqrt(tau_compare2^2 + z_compare^2);
    end
end

[rss_sort,ind] = sort(rss_compare, 'ascend');
r_compare_sort = r_compare(ind,:);
z_r_out = r_compare_sort(1,1);

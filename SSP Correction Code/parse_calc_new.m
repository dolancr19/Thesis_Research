clear variables
load calc_new
load calc_r
plot_calc = zeros(20,85);

for aa = 1:10
    ff = aa + 3*(aa-1);
    bb = 4*aa;
    gg = aa + (aa-1);
    cc = 1;
    while calc_new(ff,cc) > 0
        cc = cc+1;
    end
    cc = cc-1;
    dd = 10:10:cc;
    if dd(end) ~= cc
        dd = [dd cc];
    else
    end
        for ee = 1:length(dd)
            plot_calc(gg,ee) = abs(calc_new(bb,dd(ee))-calc_r(aa,ee));
            plot_calc(gg+1,ee) = calc_new(ff,dd(ee));
            
        end
    %save('calc_new_r', 'calc_new_r', 'cal_new_z')
end
figure
plot(plot_calc(1,1:60),plot_calc(2,1:60),'r',plot_calc(3,1:61),plot_calc(4,1:61),'g',plot_calc(5,1:61),plot_calc(6,1:61),'b',plot_calc(7,1:63),plot_calc(8,1:63),'m',plot_calc(9,1:64),plot_calc(10,1:64),'k',plot_calc(11,1:66),plot_calc(12,1:66),'r--',plot_calc(13,1:69),plot_calc(14,1:69),'g--',plot_calc(15,1:73),plot_calc(16,1:73),'b--',plot_calc(17,1:78),plot_calc(18,1:78),'m--',plot_calc(19,1:84),plot_calc(20,1:84),'k--')
title('Calculated Position Error vs Depth with Angles of Transmission between 90 and 45 Degrees')
set(gca,'YDir','Reverse')
xlabel('Horizontal Error (m)')
ylabel('Depth (m)')
legend('90', '85', '80', '75', '70', '65', '60', '55', '50', '45', 'location', 'southeast')

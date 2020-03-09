for jj = 1:length(depth_under)        
    kk = 1 + 4*(jj-1);
    ii = 1;
    while calc_new(kk,ii) < depth_under(1,jj)
        ii = ii + 1;
    end
diff1 = abs(calc_new(kk,ii)-depth_under(1,jj));
diff2 = abs(calc_new(kk,ii-1)-depth_under(1,jj));
if diff1 < diff2
    calc_new_out(jj,1) = calc_new(kk+1,ii);
    calc_new_out(jj,2) = calc_new(kk+3,ii);
else 
    calc_new_out(jj,1) = calc_new(kk+1,ii-1);
    calc_new_out(jj,2) = calc_new(kk+3,ii-1);
end
end

for ii = 1:length(calc_r(:,1))
    rnv_lat_calc_ray(ii) = wg_lat_under(ii) - (calc_r(ii,1)*sind(azi_under(ii)))/111000;
    rnv_lon_calc_ray(ii) = wg_lon_under(ii) - (calc_r(ii,1)*cosd(azi_under(ii)))/111000;
end

plot(rnv_lon_calc_ray, rnv_lat_calc_ray, 'r--', rnv_lon_under, rnv_lat_under, 'b--', wg_lon_under, wg_lat_under, 'k')

for jj = 1:length(rnv_lon_under)
    diff_lon(jj) = abs(rnv_lon_under(jj)-rnv_lon_calc(jj));
end
for kk = 59:75
    azi_under(kk) = -1*azi_under(kk);
end

figure
plot(rnv_lon_calc_ray, rnv_lat_calc_ray, 'r--o', rnv_lon_under, rnv_lat_under, 'b--o', wg_lon_under, wg_lat_under, 'ko')
xticks([-64.548:.001:-64.536]);
yticks([32.704:.001:32.714]);
title('SSP Calculated Position vs Actual Position')
xlabel('Longitude (degrees)')
ylabel('Latitude (degrees)')
legend('SSP Calculated Position', 'Actual Position', 'Wave Glider Position', 'location', 'northeast')

for kk = 1:length(hr_under)
    diff(kk,1) = kk;
    diff(kk,2) = abs(hr_under(kk)-calc_r(kk,1));
end
plot(diff(:,1), diff(:,2))

figure
plot(diff(:,1), diff(:,2))
title('Difference Between SSP Calculated Position and Actual Position')
xlabel('Sample Number')
ylabel('Error (meters)')
legend('SSP Calculated Position', 'Actual Position', 'Wave Glider Position', 'location', 'northeast')

ii = 1;
while usb_sub2.t(ii,1)-rnv.t(1,1) < 0
    ii = ii + 1;
end
ind = ii

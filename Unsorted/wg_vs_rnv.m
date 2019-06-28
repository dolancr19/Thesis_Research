wg_lat_over = interp1(wg.t, wg.lat, geopos_rnv.t);
wg_lon_over = interp1(wg.t, wg.lon, geopos_rnv.t);
diff_y = zeros(length(wg_lat_over),1);
diff_x = zeros(length(wg_lat_over),1);
hr = zeros(length(wg_lat_over),1);
angle = zeros(length(wg_lat_over),1);
for jj = 1:length(wg_lat_over)
k = 111000;
diff_y(jj) = (wg_lat_over(jj) - geopos_rnv.lat(jj))*k;
diff_x(jj) = (wg_lon_over(jj) - geopos_rnv.lon(jj))*k;
hr(jj) = sqrt((diff_x(jj))^2 + (diff_y(jj))^2);
angle(jj) = atand(diff_y(jj)/diff_x(jj));
end
save('wg_vs_geopos_rnv.mat', 'angle', 'diff_x', 'diff_y', 'wg_lat_over', 'wg_lon_over', 'hr') 

for mm = 1:round(length(hr)/300)
    nn = 1 + (300*(mm-1));
    hr_under(mm) = hr(nn);
    
end

for ii = 1:length(diff_y)
    azi(ii) = atand(diff_y(ii)/diff_x(ii));
end
save('wg_vs_rnv.mat', 'azi','ele', 'diff_x', 'diff_y', 'wg_lat_over', 'wg_lon_over', 'hr') 

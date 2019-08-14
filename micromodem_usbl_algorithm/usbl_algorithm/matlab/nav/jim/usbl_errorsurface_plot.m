% USBL error surface, treating error between actual and estimate (az,el) as
% a "Euclidean" sort of measure
error_surface_deg = NaN*ones(length(azimuth0_deg_vec), length(elevation0_deg_vec));
for azi=(1:length(azimuth0_deg_vec))
    azimuth0_deg = azimuth0_deg_vec(azi);
    for eli=(1:length(elevation0_deg_vec))
        elevation0_deg = elevation0_deg_vec(eli);
        az_hat_deg = az_hat_deg_matrix(azi,eli);
        el_hat_deg = el_hat_deg_matrix(azi,eli);
        
        az_error_deg = az_hat_deg - azimuth0_deg_vec(azi);
        el_error_deg = el_hat_deg - elevation0_deg_vec(eli);
        error_deg = sqrt(az_error_deg^2 + el_error_deg^2);
        error_surface_deg(azi,eli) = error_deg;
    end
end
hp3 = imagesc(azimuth0_deg_vec, elevation0_deg_vec, error_surface_deg);
grid on;
hx3=xlabel('azimuth (deg)');
hy3=ylabel('elevation (deg)');
if isinf(snr_dB)
	title_str3 = sprintf('USBL test (SNR=Inf, alg=%s)', search_algorithm_name);
else
	title_str3 = sprintf('USBL test (SNR=%.0f dB)', snr_dB, search_algorithm_name);
end
ht3=title(title_str3);
xlim(90*[-1 1]);
ylim(90*[-1 1]);
caxis([0 colorbar_max_error_deg]);
hc3=colorbar;
set(get(hc3,'YLabel'),'String','Error ("Euclidean" degrees)','FontSize',12);
set(hx3,'FontSize',16);
set(hy3,'FontSize',16);
set(ht3,'FontSize',24);
set(gca,'FontSize',16);
if isinf(snr_dB)
    print_command3 = sprintf('print -dpng usbl_%s_snrInf_errorsurface.png', search_algorithm_name);
else
    print_command3 = sprintf('print -dpng usbl_%s_snr%.0f_errorsurface.png', search_algorithm_name, snr_dB);
end
eval(print_command3);

% USBL error plot with little lines between actual and estimate (az,el)
for azi=(1:length(azimuth0_deg_vec))
    azimuth0_deg = azimuth0_deg_vec(azi);
    for eli=(1:length(elevation0_deg_vec))
        elevation0_deg = elevation0_deg_vec(eli);
        az_hat_deg = az_hat_deg_matrix(azi,eli);
        el_hat_deg = el_hat_deg_matrix(azi,eli);
        
        hold on;
        hp2a = plot(az_hat_deg, el_hat_deg, 'ro', azimuth0_deg, elevation0_deg, 'go');
        hp2b = line([az_hat_deg; azimuth0_deg], [el_hat_deg; elevation0_deg]);
        hold off;
    end
end
grid on;
hx2=xlabel('azimuth (deg)');
hy2=ylabel('elevation (deg)');
hl2=legend([hp2a(2) hp2a(1)], 'actual', 'estimated');
if isinf(snr_dB)
	title_str2 = sprintf('USBL test (SNR=Inf, alg=%s)', search_algorithm_name);
else
	title_str2 = sprintf('USBL test (SNR=%.0f dB)', snr_dB, search_algorithm_name);
end
ht2=title(title_str2);
xlim(90*[-1 1]);
ylim(90*[-1 1]);
set(hx2,'FontSize',16);
set(hy2,'FontSize',16);
set(ht2,'FontSize',24);
set(gca,'FontSize',16);
if isinf(snr_dB)
    print_command2 = sprintf('print -dpng usbl_%s_snrInf_errorlines.png', search_algorithm_name);
else
    print_command2 = sprintf('print -dpng usbl_%s_snr%.0f_errorlines.png', search_algorithm_name, snr_dB);
end
eval(print_command2);

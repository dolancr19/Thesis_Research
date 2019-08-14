% usbl_test_jwp.m     script for testing usbl solutions

path('../../misc',path);

plot_shifted_ires=0;
search_algorithm_name = 'bruteforce';

% actual azimuth,elevation for phase/delay-shifting signals
azimuth0_deg_vec   =  (-85:5:85);
elevation0_deg_vec =  (-85:5:85);
% estimated azimuth,elevation from USBL algorithm
az_hat_deg_matrix = NaN*ones(length(azimuth0_deg_vec), length(elevation0_deg_vec));
el_hat_deg_matrix = NaN*ones(length(azimuth0_deg_vec), length(elevation0_deg_vec));


% FIXME: I'd like these to be the search vectors, but instead it's done in
% terms of axvec and ayvec, which are components of the direction/normal
% vector.
%azimuth_search_deg_vec   = (-(90-1):1:(90-1));
%elevation_search_deg_vec = (-(90-1):1:(90-1));


fc_Hz = 15000; % carrier
fs_Hz = 80000 ; % passband
bw_Hz = 5000; 
nnew = 2 ;
fb_Hz = bw_Hz*nnew;
snr_dB = Inf;
c_ms = 1500;
dthresh = 10 ;
pown = 20 ;

lambda_m = c_ms/fc_Hz;
k0_1m = 2*pi/lambda_m;

d_m = (lambda_m/2)/2;
% square array
%xyz_m  = (2*d_m)*[[0 0 0]' [1 0 0]' [1 1 0]' [0 1 0]'];
% diagonal array so that azimuth and elevation are easily separated
xyz_m  = d_m*[[0 1 0]' [1 0 0]' [0 -1 0]' [-1 0 0]'];

f1_Hz = fc_Hz - bw_Hz/2;
f2_Hz = fc_Hz + bw_Hz/2;
T0_s = 10e-3;
y0 = lfm0(f1_Hz, f2_Hz, T0_s, fs_Hz, 0).';
iv0 = p2b(y0,[fb_Hz fc_Hz]/fs_Hz) ;
z0 = zeros(1e3,1);
iv1 = [z0;iv0;z0];
h_mfd = conj(iv0); % flipud() is done in mfd1()

for azi=(1:length(azimuth0_deg_vec))
    azimuth0_deg = azimuth0_deg_vec(azi);
    azimuth0_rad = (pi/180)*azimuth0_deg;
    fprintf('azi=%3d of %d\n', azi, length(azimuth0_deg_vec));    
    for eli=(1:length(elevation0_deg_vec))
        elevation0_deg = elevation0_deg_vec(eli);
        elevation0_rad = (pi/180)*elevation0_deg;

        % normal vector is FROM source TO receiver, hence minus sign
        normal_vector = -[tan(azimuth0_rad) tan(elevation0_rad) 1]';
        unit_normal_vector = normal_vector/sqrt(normal_vector'*normal_vector);
        k_vec_1m = k0_1m * unit_normal_vector;

        iv = delay_and_phase_shift_array_signal(xyz_m, fc_Hz, unit_normal_vector, iv1, c_ms, fb_Hz, snr_dB);

        [ire,ps,dv,p] = mfd1(iv,h_mfd,dthresh,pown);
        [peak_value,peak_sample] = max(abs(ire(:,1)));

        X = transpose(ire(peak_sample,:));

        if plot_shifted_ires
            figure(1); clf;
            aX = angle(X);
            az_hat0_deg = (180/pi)*asin(-(aX(4)-aX(2))/(k0_1m*2*d_m));
            el_hat0_deg = (180/pi)*asin(-(aX(3)-aX(1))/(k0_1m*2*d_m));
            title_str1 = sprintf('az=%.1f deg, azhat0=%.1f deg; el=%.1f deg, elhat0=%.1f deg', azimuth0_deg, az_hat0_deg, elevation0_deg, el_hat0_deg);
            xlim_vec=peak_sample+[-25 25];
            ylim_vec=8*[-1 1];
            h1_vec = usbl_ire_plots(ire, peak_sample, title_str1, xlim_vec, ylim_vec);
            %pause;
        end; % plot_shifted_ires
        
        ax_search_vec = (-1:0.01:+1)';
        ay_search_vec = (-1:0.01:+1)';
        [az_hat_deg, el_hat_deg, ximax, yimax] = usbl_search_brute_force(X, xyz_m, k0_1m, ax_search_vec, ay_search_vec);
        az_hat_deg_matrix(azi,eli) = az_hat_deg;
        el_hat_deg_matrix(azi,eli) = el_hat_deg;
    end
end
if isinf(snr_dB)
    usbl_data_filename = sprintf('usbl_%s_snrInf_data.mat', search_algorithm_name);
else
    usbl_data_filename = sprintf('save usbl_%s_snr%.0f_data.mat', search_algorithm_name, snr_dB);
end
save(usbl_data_filename);

%load(usbl_data_filename);

% USBL error plot with little lines between actual and estimate (az,el)
figure(2); clf
usbl_errorlines_plot;

% USBL error surface, treating error between actual and estimate (az,el) as
% a "Euclidean" sort of measure
figure(3); clf
colorbar_max_error_deg=8;
usbl_errorsurface_plot;

return;


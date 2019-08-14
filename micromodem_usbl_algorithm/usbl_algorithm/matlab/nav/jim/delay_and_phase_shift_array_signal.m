% function to create array data for a 3-d array
% with arbitrary element spacing.
%
% Inputs:
%       xyz_m            3-d vector of sensor locations:
%                            [x1 y1 z1; x2 y2 z2; ... ; xn yn zn]
%       fc_Hz            signal's center frequency (Hz)
%       unit_normal_vector  [dx dy dz]' unit normal vector of plane waves on array
%       y0_bb            undelayed baseband signal
%       c_ms             sound speed (m/s)
%       fs_Hz            sample frequency (Hz)
%       snr_dB           signal-to-noise ratio (dB)
%
% function y = twodasig(xy_m,fc_Hz,unit_normal_vector,y0,c_ms,fs_Hz,snr_dB)

function y = twodasig(xyz_m,fc_Hz,unit_normal_vector,y0_bb,c_ms,fb_Hz,snr_dB)

numsensors = size(xyz_m,2);
Tb_s = 1/fb_Hz;
sigtime_s = length(y0_bb)*Tb_s;
numpoints0 = floor(sigtime_s/Tb_s);
pad_length = 100;
numpoints1 = 2*pad_length + numpoints0;

y = NaN*ones(numpoints1, numsensors);

% sensor 1 is the reference pathlength and phase
xyz_m_origin = xyz_m(:,1);
y_ii = zeros(numpoints1,1);
n0 = pad_length + 0;
n1 = n0 + numpoints0 - 1;
y_ii(n0:n1) = y0_bb;
y(:,1) = y_ii;

% sensors 2 through N
for ii=(2:numsensors)
    delta_pathlength_m = unit_normal_vector'*(xyz_m(:,ii)-xyz_m_origin);
    tau0_s = delta_pathlength_m/c_ms;
    delta_N = fix(tau0_s/Tb_s);
    tau1_s = tau0_s - delta_N*Tb_s;
    phase_shift_rad = tau1_s * (2*pi*fc_Hz);
    %fprintf('ii=%d, ds_m=%f m, tau0=%f ms, tau1=%f ms, delta_N=%d, phase_shift %f deg\n',ii,delta_pathlength_m,1e3*tau0_s,1e3*tau1_s,delta_N,(180/pi)*phase_shift_rad);
    
    y_ii = zeros(numpoints1,1);
    n0 = pad_length + delta_N;
    n1 = n0 + numpoints0 - 1;
    y_ii(n0:n1) = y0_bb*exp(sqrt(-1)*phase_shift_rad);
    y(:,ii) = y_ii;
end

% noise
if isinf(snr_dB)
    return;
end
noise_mag = 10^(-snr_dB/20);
for ii=(1:numsensors)
    noise0_ii = (noise_mag/sqrt(2)) * (randn(numpoints1, 1) + sqrt(-1)*randn(numpoints1, 1));
    noise_ii  = noise0_ii - mean(noise0_ii);
    y(:,ii) = y(:,ii) + noise_ii;
end

return;

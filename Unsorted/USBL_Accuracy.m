Gamma = 0;%slant range
alpha = 0;%azimuth
gamma = 0;%elevation

Phi = 0;%ship roll
theta = 0;%ship pitch
psi = 0;%ship heading

Phi_prime = 0;%transducer offet in roll
theta_prime = 0;%transducer offset in pitch
psi_prime = 0;%transducer offset in heading

%pt = position relative to USBL transducer in the reference frame of the
%transducer
%pt(Gamma,alpha,gamma);

%Rtb = rotation matrix from transducer reference frame to ship reference frame
%Rtb(Phi_prime,theta_prime,psi_prime);

%Rbn = rotation matrix from ship reference frame to North-East-Down
%reference frame
%Rbn(Phi,theta,psi);

%pn = position of AUV measured relative to the USBL transducer
%pn = Rbn(Phi,theta,psi)*Rtb(Phi_prime,theta_prime,psi_prime)*pt(Gamma,alpha,gamma);

clear variables

x = 1000; %m
y = 1000; %m
z = [50, 100, 500, 1000, 3000]; %m

inc_error = .12; %degrees at 20 dB
azi_error = .12; %degrees at 20 dB
GPS_error = .2; %m

for ii = 1:length(z)
    r = sqrt(x^2+y^2+z(ii)^2)
    inc = acosd(z(ii)/r);
    azi = atand(y/x);
    x_plus = (r*sind(inc - inc_error/2)*cosd(azi - azi_error/2));
    x_minus = (r*sind(inc + inc_error/2)*cosd(azi + azi_error/2));
    delta_x = x_plus - x_minus + GPS_error/2;

    y_plus = (r*sind(inc + inc_error/2)*sind(azi + azi_error/2));
    y_minus = (r*sind(inc - inc_error/2)*sind(azi - azi_error/2));
    delta_y = y_plus - y_minus + GPS_error/2;

    %pu = position uncertainty of the AUV
    pu = sqrt(delta_x^2 + delta_y^2);
    
    out(ii,:) = [ii pu]
end

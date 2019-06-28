clear all

x = 50; %m
y = 50; %m
z = 50; %m

inc_error = 0; %degrees at 20 dB
azi_error = .12; %degrees at 20 dB
GPS_error = 0; %m

Razi_plus = rotz(-azi_error/2);
Razi_minus = rotz(azi_error/2);
pos = [x y z];
azi_prime_plus = Razi_plus * pos.'
azi_prime_minus = Razi_minus * pos.'
azi_x_prime = abs(azi_prime_plus(1) - azi_prime_minus(1))
azi_y_prime = abs(azi_prime_plus(2) - azi_prime_minus(2))

Rinc_plus = roty(inc_error/2);
Rinc_minus = roty(-inc_error/2);
inc_prime_plus = Rinc_plus * pos.'
inc_prime_minus = Rinc_minus * pos.'
inc_x_prime = abs(inc_prime_plus(1) - inc_prime_minus(1))
inc_y_prime = abs(inc_prime_plus(2) - inc_prime_minus(2))

horiz_error = sqrt(azi_x_prime^2 + azi_y_prime^2)
for ii = 1:length(z)
    r = sqrt(x^2+y^2+z(ii)^2);
    inc = acosd(z(ii)/r);
    azi = atand(y/x);
    x_plus = (r*sind(inc - inc_error/2)*cosd(azi - azi_error/2));
    x_minus = (r*sind(inc + inc_error/2)*cosd(azi + azi_error/2));
    delta_x = x_plus - x_minus + GPS_error/2

    y_plus = (r*sind(inc + inc_error/2)*sind(azi + azi_error/2));
    y_minus = (r*sind(inc - inc_error/2)*sind(azi - azi_error/2));
    delta_y = y_plus - y_minus + GPS_error/2

    %pu = position uncertainty of the AUV
    pu = sqrt(delta_x^2 + delta_y^2);
    
    out(ii,:) = [ii pu]
end

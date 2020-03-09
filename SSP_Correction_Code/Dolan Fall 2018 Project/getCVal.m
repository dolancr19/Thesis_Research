%% Helper Function: getCVal;

function [c,cz,czz]=getCVal(r,z,cGridVal)
%return a c value in m/s for a given r,z value:
%search grid
c=0; % soundspeed
cz=0; % slope
czz=0; % second derivative

% NOTE: You can pass in some structure, cGridVal, to make this function three
% lines, or you can do an interpolation here. What we strongly suggest you
% DO NOT do is find the derivative and interpolation object every single
% time because that is a redundant amount of calculations.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% INSERT YOUR CODE HERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% TO GET soundspeed and soundspeed gradient for r,z  %%%%%%%%%%%%

% cVal is range-independent, don't need r;
c = cGridVal(z);

dz = 0.5;

if z>dz
    c_minus = cGridVal(z-dz);
else
    c_minus = c;
end
c_plus = cGridVal(z+dz);

cz = (c - c_minus)./dz;
cz_plus = (c_plus-c)./dz;

czz = (cz_plus - cz)./dz;

if ~isreal(z)
    warning('Depth is imaginary:')
    display(z)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
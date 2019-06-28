%% Helper Function: getDepth

function [D,ang]=getDepth(r,dGridVal,dzVal)

D=5000;
ang=0;
% r- range value in m
% dGridVal- cell array of depth griddedInterpolants

% NOTE: You can pass in some structure, dGridVal, to make this function one
% line, or you can do an interpolation here. What we strongly suggest you
% DO NOT do is find the derivative and interpolation object every single
% time because that is a redundant amount of calculations.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% INSERT YOUR CODE HERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% TO GET DEPTH D and bottom angle ang at cur r, dVal  %%%%%%%%%%

D = dGridVal(r);
dzVal = (D(2:end)-D(1:end-1))/(r(2:end)-r(1:end-1));
ang = atan(dzVal);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
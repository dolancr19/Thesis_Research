function [x,y,z] = vrot(X,Y,Z,varargin)
%[x,y,z] = vrot(X,Y,Z,roll,pitch,heading,{'in'|'out'})
%  'in' this steps in an axis: yaw the original axis about Z, pitch about Y, roll about X.
%  'out' this steps out an axis: -roll about x, -pitch about y, -yaw about z
%
%   [x,y,z]' = R('in')[X,Y,Z]'
%   [X,Y,Z]' = R('out')[x,y,z]'
%
%   where R('in') = R('out').'
%
%[x,y,z] = vrot(X,Y,Z,R,{'in'|'out'})  % specify the rotation matrix

% Revision History
% 2004-?-?      mvj    1.0    created
% 2005-03-16    mvj    1.0    added option to specify R   
% 2005-03-30    mvj    1.0    matrix x,y,z arguments supported, assumes row-wise correspondence


% parse inputs
if nargin == 7
  roll = varargin{1};
  pitch = varargin{2}; 
  heading = varargin{3}; 
  in_or_out = varargin{4};
    
  % precalculate trigonometric quantities
  cf = cos(roll(:));
  sf = sin(roll(:));
  ct = cos(pitch(:));
  st = sin(pitch(:));
  cy = cos(heading(:));
  sy = sin(heading(:));
  
  % vector versions of each element of the rotation matrix for each rph triple
  R11 = cy.*ct;
  R12 = -sy.*cf+cy.*st.*sf;
  R13= sy.*sf+cy.*cf.*st;
  R21 = sy.*ct;
  R22 = cy.*cf+sf.*st.*sy;
  R23 = -cy.*sf+st.*sy.*cf;
  R31 = -st;
  R32 = ct.*sf;
  R33 = ct.*cf;

elseif nargin == 5
  R = varargin{1};
  in_or_out = varargin{2};
  
  R11 = R(1,1);
  R12 = R(1,2);
  R13 = R(1,3);
  R21 = R(2,1);
  R22 = R(2,2); 
  R23 = R(2,3);
  R31 = R(3,1);
  R32 = R(3,2);
  R33 = R(3,3);

end

if strcmp(in_or_out,'out')  % [X,Y,Z]' = R('out')[x,y,z]'  % this is what rot_matrix returns (-roll abot x0, -pitch
                                                           % about y1, -yaw about z2)

  % do nothing
  
elseif strcmp(in_or_out,'in') % [x,y,z]' = R('in')[X,Y,Z]' % yaw about Z3, pitch about y2, roll about x1 

  % transpose
  [R11,R12,R13,R21,R22,R23,R31,R32,R33] = deal(R11,R21,R31,R12,R22,R32,R13,R23,R33);
  
else
  error('specify({''in''|''out''})');
end

% now do the rotation in vector form
[M,N] = size(X);
if N > 1 & length(R11) > 1
  R11 = repmat(R11(:),N,1);
  R12 = repmat(R12(:),N,1);
  R13 = repmat(R13(:),N,1);
  R21 = repmat(R21(:),N,1);
  R22 = repmat(R22(:),N,1);
  R23 = repmat(R23(:),N,1);
  R31 = repmat(R31(:),N,1);
  R32 = repmat(R32(:),N,1);
  R33 = repmat(R33(:),N,1);
  x = R11.*X(:)+R12.*Y(:)+R13.*Z(:);
  y = R21.*X(:)+R22.*Y(:)+R23.*Z(:);
  z = R31.*X(:)+R32.*Y(:)+R33.*Z(:);
  x = reshape(x,M,N);
  y = reshape(y,M,N);
  z = reshape(z,M,N);
else
  x = R11.*X+R12.*Y+R13.*Z;
  y = R21.*X+R22.*Y+R23.*Z;
  z = R31.*X+R32.*Y+R33.*Z;
end 

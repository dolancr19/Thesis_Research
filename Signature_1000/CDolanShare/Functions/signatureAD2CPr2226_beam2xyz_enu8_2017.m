function [ Data, Config, T_beam2xyz ] = signatureAD2CPr2226_beam2xyz_enu8_2017( Data, Config, mode, twoZs )

% updated 8/2014 for latest Nortek fireware and software

% ** This function has been altered to allow proper rotation.  The Config
% structure has incorrect theta values for the alignment of the beams with
% the instrument 'X'.  The alteration (Line 52) resets the
% beamVectorsSpherical matrix such that ThetaBeam1 = 0, ThetaBeam2 = -90,
% ThetaBeam3 = 180, ThetaBeam4 = 90. 8/12/2014

% ** Additional modification was applied to the ENU transformation.  16.1 
% degrees was added to the heading to account for the difference in true 
% magnetic north as of 8/15/2014. 
% signXYZ, which determines if the ADCP is downward or upward facing based 
% on the status bits did not produce the correct transformation.  Therefore,
% flipping the E and Up axis was preformed after the ENU transformation as
% opposed to before to the YZ axes.



if nargin == 3
	twoZs = 0;
end

ad2cpstr = '';

if strcmpi( mode, 'avg' )
	dataModeWord = 'Average';
	configModeWord = 'avg';
elseif strcmpi( mode, 'burst' )
	dataModeWord = 'Burst';
   BTstr = '';
	configModeWord = 'burst';
elseif strcmpi( mode, 'bt' )
	dataModeWord = 'BurstBT';
   BTstr = 'BT';
	configModeWord = 'bt';
   Config.bt_coordSystem = 'beam';  % Add field for this version of the data
   Config.bt_nCells = 1;            % Correct invalid value for nCells
   Data.BurstBT_TimeStamp = Data.IBurst_MatlabTimeStamp(1:length(Data.BurstBT_VelBeam1)); % Fix missing time stamp
end

% make the assumption the beam mapping is the same for all measurements in a data file
activeBeams = Data.( [dataModeWord '_Physicalbeam'] )( 1, : );
activeBeams = activeBeams(find(activeBeams > 0));
numberOfBeams = length( activeBeams );
if numberOfBeams <= 2
	print( 'Transformations require at least 3 active beams.' )
	T_beam2xyz = nan;
	return
end
% assume max number of beams involved is 4, extra rows removed later
beamVectorsSpherical = zeros( 4, 3 );

for i = activeBeams
	beamVectorsSpherical( i, : ) = [ 1, ...
		Config.( [ad2cpstr 'beamConfiguration' num2str( i ) '_theta' ] ), ...
		Config.( [ad2cpstr 'beamConfiguration' num2str( i ) '_phi' ] ) ];
end


% beam locations based on true theta values.  positive x points in the 
% direction of the boat
beamVectorsSpherical=[1 25 0; 1 25 -90; 1 25 180; 1 25 90];

disp(beamVectorsSpherical)

if numberOfBeams == 3
	% for a three beam system, translate the beam vectors expressed in spherical coordinates
	% into beam vectors in Cartesian coordinates

	% first transform from spherical to cartesian coordinates
	for i = activeBeams
		beamVectorsCartesian( i, : ) = [ ...
			sind( beamVectorsSpherical( i, 2 ) ) * cosd( beamVectorsSpherical( i, 3 ) ), ...
			sind( beamVectorsSpherical( i, 2 ) ) * sind( beamVectorsSpherical( i, 3 ) ), ...
			cosd( beamVectorsSpherical( i, 2 ) ) ];
    end


	cartesianTransformCheck = sum( beamVectorsCartesian.^2, 2 );

	% remove any extra rows for inactive beams
	beamVectorsCartesian( cartesianTransformCheck == 0, : ) = [];

    
	T_beam2xyz = inv( beamVectorsCartesian );
  
elseif numberOfBeams == 4
	if twoZs == 0
		% for a four beam system, translate the beam vectors expressed in spherical coordinates
		% into beam vectors in Cartesian coordinates, using only three basis vectors
		for i = 1:numberOfBeams
			beamVectorsCartesian( i, : ) = [ ...
				sind( beamVectorsSpherical( i, 2 ) ) * cosd( beamVectorsSpherical( i, 3 ) ), ...
				sind( beamVectorsSpherical( i, 2 ) ) * sind( beamVectorsSpherical( i, 3 ) ), ...
				cosd( beamVectorsSpherical( i, 2 ) ) ];
		end
		cartesianTransformCheck = sum( beamVectorsCartesian.^2, 2 );

		% pseudo inverse needs to be used because beamVectorsCartesian isn't square
		T_beam2xyz = pinv( beamVectorsCartesian );

	else
		% this section makes two estimates of the vertical velocity
		for i = 1:numberOfBeams
			if i == 1 | i == 3
				beamVectorsCartesianzz( i, : ) = [ ...
					beamVectorsSpherical( i, 1 ) * sind( beamVectorsSpherical( i, 2 ) ) * cosd( beamVectorsSpherical( i, 3 ) ), ...
					-1 * beamVectorsSpherical( i, 1 ) * sind( beamVectorsSpherical( i, 2 ) ) * sind( beamVectorsSpherical( i, 3 ) ), ...
					beamVectorsSpherical( i, 1 ) * cosd( beamVectorsSpherical( i, 2 ) ), ...
					0 ];
			else
				beamVectorsCartesianzz( i, : ) = [ ...
					beamVectorsSpherical( i, 1 ) * sind( beamVectorsSpherical( i, 2 ) ) * cosd( beamVectorsSpherical( i, 3 ) ), ...
					-1 * beamVectorsSpherical( i, 1 ) * sind( beamVectorsSpherical( i, 2 ) ) * sind( beamVectorsSpherical( i, 3 ) ), ...
					0, ...
					beamVectorsSpherical( i, 1 ) * cosd( beamVectorsSpherical( i, 2 ) ) ];
            end

		end
		cartesianTransformCheck = sum( beamVectorsCartesianzz.^2, 2 );

		% there should be an inverse for this, no pseudoinverse needed
		T_beam2xyz = inv( beamVectorsCartesianzz );


        
		% Can also add in a row for the error velocity calculation, 
		% as the difference between the two vertical velcoities
		% T_beam2xyz( end + 1, : ) = T_beam2xyzz( 3, : ) - T_beam2xyzz( 4, : );
		% note the addition of the error velocity row changes the inversion of the matrix, it
		% needs to be removed to recover the xyz2beam matrix.
    end
    

end

% verify we're not already in 'xyz'
if strcmpi( Config.( [ ad2cpstr configModeWord '_coordSystem' ] ), 'xyz' )
	disp( 'Velocity data is already in xyz coordinate system.' )
	return
end

xAllCells = zeros( length( Data.( [ dataModeWord '_TimeStamp' ] ) ), Config.( [ad2cpstr  configModeWord '_nCells' ] ) );
yAllCells = zeros( length( Data.( [ dataModeWord '_TimeStamp' ] ) ), Config.( [ad2cpstr  configModeWord '_nCells' ] ) );
zAllCells = zeros( length( Data.( [ dataModeWord '_TimeStamp' ] ) ), Config.( [ad2cpstr  configModeWord '_nCells' ] ) );
if twoZs == 1
	z2AllCells = zeros( length( Data.( [ dataModeWord '_TimeStamp' ] ) ), Config.( [ad2cpstr  configModeWord '_nCells' ] ) );
end

xyz = zeros( size( T_beam2xyz, 2 ), length( Data.( [ dataModeWord '_TimeStamp' ] ) ) );
beam = zeros( size( T_beam2xyz, 2 ), length( Data.( [ dataModeWord '_TimeStamp' ] ) ) );
for nCell = 1:Config.( [ad2cpstr  configModeWord '_nCells' ] )
	for i = 1:numberOfBeams
		beam( i, : ) = Data.( [ dataModeWord '_' 'VelBeam' num2str( Data.( [ dataModeWord '_Physicalbeam' ] )( 1, i ) ) ] )( :, nCell )';
	end
	xyz = T_beam2xyz * beam;
	xAllCells( :, nCell ) = xyz( 1, : )';	
	yAllCells( :, nCell ) = xyz( 2, : )';
	zAllCells( :, nCell ) = xyz( 3, : )';
	if twoZs == 1
		z2AllCells( :, nCell ) = xyz( 4, : )';
	end
end

Config.( [ad2cpstr   configModeWord '_coordSystem' ] ) = 'xyz';
Data.( [ dataModeWord '_VelX' ] ) = xAllCells;
Data.( [ dataModeWord '_VelY' ] ) = yAllCells;

if twoZs == 1
	Data.( [ dataModeWord '_VelZ1' ] ) = zAllCells;
	Data.( [ dataModeWord '_VelZ2' ] ) = z2AllCells;
else
	Data.( [ dataModeWord '_VelZ' ] ) = zAllCells;
end





% verify we're not already in 'enu'
if strcmpi( Config.( [ad2cpstr   configModeWord '_coordSystem' ] ), 'enu' )
	disp( 'Velocity data is already in enu coordinate system.' )
	return
end

K = 3;
EAllCells = zeros( length( Data.( [dataModeWord  '_TimeStamp' ] ) ), Config.( [ad2cpstr   configModeWord '_nCells' ] ) );
NAllCells = zeros( length( Data.( [dataModeWord  '_TimeStamp' ] ) ), Config.( [ad2cpstr   configModeWord '_nCells' ] ) );
UAllCells = zeros( length( Data.( [dataModeWord  '_TimeStamp' ] ) ), Config.( [ad2cpstr   configModeWord '_nCells' ] ) );
if twoZs == 1
	U2AllCells = zeros( length( Data.( [dataModeWord  '_TimeStamp' ] ) ), Config.( [ad2cpstr  configModeWord '_nCells' ] ) );
   K = 4;
end

Name = ['X','Y','Z'];
ENU = zeros( K, Config.([ad2cpstr   configModeWord '_nCells' ]));
xyz = zeros( K, Config.([ad2cpstr   configModeWord '_nCells' ]));
for sampleIndex = 1:length(Data.( [dataModeWord  '_Error' ]));
   if (bitand(bitshift(uint32(Data.( [dataModeWord  '_Status' ])(sampleIndex)), -25),7) == 5)
      signXYZ=[1 -1 -1 -1];
   else
      signXYZ=[1 1 1 1];
   end
   
   signXYZ = [1 1 1 1]; % -- changed: assume upward facing ADCP and account for
   % downward facing coordinate system after ENU transformation
   
   hh = pi*(Data.([dataModeWord  '_Heading'])(sampleIndex)+16.1-90)/180; % originally just -90
   % add 180 to see what way the boat is going (bt away from boat direction
   % 16.1 is magnetic vs true north as of 8/15/2014
   pp = pi*Data.([dataModeWord  '_Pitch'])(sampleIndex)/180;
   rr = pi*Data.([dataModeWord  '_Roll'])(sampleIndex)/180;

   % Make heading matrix
   H = [cos(hh) sin(hh) 0; -sin(hh) cos(hh) 0; 0 0 1];

   % Make tilt matrix
   P = [cos(pp) -sin(pp)*sin(rr) -cos(rr)*sin(pp);...
         0             cos(rr)          -sin(rr);  ...
         sin(pp) sin(rr)*cos(pp)  cos(pp)*cos(rr)];

   % Make resulting transformation matrix
   xyz2enu = H*P; 
   if (twoZs == 1)
      xyz2enu(1,3) = xyz2enu(1,3)/2;
      xyz2enu(1,4) = xyz2enu(1,3);
      xyz2enu(2,3) = xyz2enu(2,3)/2;
      xyz2enu(2,4) = xyz2enu(2,3);
      
      xyz2enu(4,:) = xyz2enu(3,:);
      xyz2enu(3,4) = 0;
      xyz2enu(4,4) = xyz2enu(3,3);
      xyz2enu(4,3) = 0;
   end

   for i = 1:K
      if (twoZs == 1) && (i >= 3)
         axs = [ Name(3) num2str((i-2),1) ];
      else
         axs = Name(i);
      end
      xyz( i, : ) = signXYZ(i) * Data.( [ dataModeWord '_Vel' axs] )( sampleIndex, : )';
   end
   ENU = xyz2enu * xyz;
   EAllCells( sampleIndex, : ) = -ENU( 1, : )';	% - downwardfacing ADCP
   NAllCells( sampleIndex, : ) = ENU( 2, : )';
   UAllCells( sampleIndex, : ) = -ENU( 3, : )'; % - downwardfacing ADCP
      if twoZs == 1
      U2AllCells( sampleIndex, : ) = -ENU( 4, : )'; % - downwardfacing ADCP
      end
%    EAllCells( sampleIndex, : ) = ENU( 1, : )';	% - removed WK
%    NAllCells( sampleIndex, : ) = ENU( 2, : )';
%    UAllCells( sampleIndex, : ) = ENU( 3, : )'; % - removed WK
%       if twoZs == 1
%       U2AllCells( sampleIndex, : ) = ENU( 4, : )'; % - removed WK
%       end      
end
Config.( [ad2cpstr   configModeWord '_coordSystem' ] ) = 'enu';
Data.( [ dataModeWord '_VelEast' ] ) = EAllCells;
Data.( [ dataModeWord '_VelNorth' ] ) = NAllCells;
if twoZs == 1
	Data.( [ dataModeWord '_VelUp1' ] ) = UAllCells;
	Data.( [ dataModeWord '_VelUp2' ] ) = U2AllCells;
else
	Data.( [ dataModeWord '_VelUp' ] ) = UAllCells;
end




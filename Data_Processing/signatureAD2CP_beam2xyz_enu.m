% Function provided by Nortek

% Modified due to modification in  Config struct 
% Line AD2CP_ was removed
% MGP, April 29th, 2015
% Fixed bug for Two VelZ solution (removed -1)

% MIT License
% 
% Copyright (c) 2017 mguerrap
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated 
% documentation files (the "Software"), to deal in the Software without restriction, including without limitation 
% the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to 
% permit persons to whom the Software is furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all copies or substantial portions of the 
% Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE 
% WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR 
% COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR 
% OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

function [ Data, Config, T_beam2xyz ] = signatureAD2CP_beam2xyz_enu( Data, Config, dataModeWord, planModeWord, twoZs, declination_correction )

% make the assumption the beam mapping is the same for all measurements in a data file
%activeBeams = Data.( [dataModeWord '_Physicalbeam'] )( 1, : );
%activeBeams = activeBeams(find(activeBeams > 0));
%numberOfBeams = length( activeBeams );
%numberOfBeams = Data.( [dataModeWord '_NBeams'] )( 1, 1 );
numberOfBeams = Data.Burst_NBeams(1,1);
if numberOfBeams <= 2
	print( 'Transformations require at least 3 active beams.' )
	T_beam2xyz = nan;
	return
end
% assume max number of beams involved is 4, extra rows removed later
beamVectorsSpherical = zeros( 4, 3 );
    
for ii = 1:numberOfBeams
	beamVectorsSpherical( ii, : ) = [ 1, ...
		Config.( ['BeamCfg' num2str( ii ) '_theta' ] ), ...
		Config.( ['BeamCfg' num2str( ii ) '_phi' ] ) ];
end

disp(beamVectorsSpherical)

if numberOfBeams == 3
	% for a three beam system, translate the beam vectors expressed in spherical coordinates
	% into beam vectors in Cartesian coordinates

	% first transform from spherical to cartesian coordinates
	for ii = activeBeams
		beamVectorsCartesian( ii, : ) = [ ...
			sind( beamVectorsSpherical( ii, 2 ) ) * cosd( beamVectorsSpherical( ii, 3 ) ), ...
			sind( beamVectorsSpherical( ii, 2 ) ) * sind( beamVectorsSpherical( ii, 3 ) ), ...
			cosd( beamVectorsSpherical( ii, 2 ) ) ];
	end

	cartesianTransformCheck = sum( beamVectorsCartesian.^2, 2 );
	% remove any extra rows for inactive beams
	beamVectorsCartesian( cartesianTransformCheck == 0, : ) = [];
	
	T_beam2xyz = inv( beamVectorsCartesian );

elseif numberOfBeams == 4
	if twoZs == 0
		% for a four beam system, translate the beam vectors expressed in spherical coordinates
		% into beam vectors in Cartesian coordinates, using only three basis vectors
		for ii = 1:numberOfBeams
			beamVectorsCartesian( ii, : ) = [ ...
				sind( beamVectorsSpherical( ii, 2 ) ) * cosd( beamVectorsSpherical( ii, 3 ) ), ...
				sind( beamVectorsSpherical( ii, 2 ) ) * sind( beamVectorsSpherical( ii, 3 ) ), ...
				cosd( beamVectorsSpherical( ii, 2 ) ) ];
		end
		cartesianTransformCheck = sum( beamVectorsCartesian.^2, 2 );

		% pseudo inverse needs to be used because beamVectorsCartesian isn't square
		T_beam2xyz = pinv( beamVectorsCartesian );

	else
		% this section makes two estimates of the vertical velocity
		for ii = 1:numberOfBeams
			if ii == 1 || ii == 3
				beamVectorsCartesianzz( ii, : ) = [ ...
					beamVectorsSpherical( ii, 1 ) * sind( beamVectorsSpherical( ii, 2 ) ) * cosd( beamVectorsSpherical( ii, 3 ) ), ...
					1 * beamVectorsSpherical( ii, 1 ) * sind( beamVectorsSpherical( ii, 2 ) ) * sind( beamVectorsSpherical( ii, 3 ) ), ...
					beamVectorsSpherical( ii, 1 ) * cosd( beamVectorsSpherical( ii, 2 ) ), ...
					0 ];
			else
				beamVectorsCartesianzz( ii, : ) = [ ...
					beamVectorsSpherical( ii, 1 ) * sind( beamVectorsSpherical( ii, 2 ) ) * cosd( beamVectorsSpherical( ii, 3 ) ), ...
					1 * beamVectorsSpherical( ii, 1 ) * sind( beamVectorsSpherical( ii, 2 ) ) * sind( beamVectorsSpherical( ii, 3 ) ), ...
					0, ...
					beamVectorsSpherical( ii, 1 ) * cosd( beamVectorsSpherical( ii, 2 ) ) ];
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
if strcmpi(Config.Burst_CoordSystem, 'xyz' )
	disp( 'Velocity data is already in xyz coordinate system.' )
	return
end

xAllCells = zeros( length( Data.( [ planModeWord '_Time' ] ) ), Config.Burst_NCells);
yAllCells = zeros( length( Data.( [ planModeWord '_Time' ] ) ), Config.Burst_NCells);
zAllCells = zeros( length( Data.( [ planModeWord '_Time' ] ) ), Config.Burst_NCells);
if twoZs == 1
	z2AllCells = zeros( length( Data.( [ planModeWord '_Time' ] ) ), Config.Burst_NCells);
end

%xyz = zeros( size( T_beam2xyz, 2 ), length( Data.( [ planModeWord '_Time' ] ) ) );
beam = zeros( size( T_beam2xyz, 2 ), length( Data.( [ planModeWord '_Time' ] ) ) );
for nCell = 1:Config.Burst_NCells
	for ii = 1:numberOfBeams
		beam( ii, : ) = Data.( [ dataModeWord '_VelBeam' num2str( ii ) ] )( :, nCell )';
	end
	xyz = T_beam2xyz * beam;
	xAllCells( :, nCell ) = xyz( 1, : )';	
	yAllCells( :, nCell ) = xyz( 2, : )';
	zAllCells( :, nCell ) = xyz( 3, : )';
	if twoZs == 1
		z2AllCells( :, nCell ) = xyz( 4, : )';
	end
end

%Config.( [ad2cpstr   dataModeWord '_CoordSystem' ] ) = 'xyz';
Data.( [ dataModeWord '_VelX' ] ) = xAllCells;
Data.( [ dataModeWord '_VelY' ] ) = yAllCells;

if twoZs == 1
	Data.( [ dataModeWord '_VelZ1' ] ) = zAllCells;
	Data.( [ dataModeWord '_VelZ2' ] ) = z2AllCells;
else
	Data.( [ dataModeWord '_VelZ' ] ) = zAllCells;
end

% verify we're not already in 'enu'
if strcmpi( Config.Burst_CoordSystem, 'enu' )
	disp( 'Velocity data is already in enu coordinate system.' )
	return
end

K = 3;
EAllCells = zeros( length( Data.( [planModeWord  '_Time' ] ) ), Config.Burst_NCells);
NAllCells = zeros( length( Data.( [planModeWord  '_Time' ] ) ), Config.Burst_NCells);
UAllCells = zeros( length( Data.( [planModeWord  '_Time' ] ) ), Config.Burst_NCells);
if twoZs == 1
	U2AllCells = zeros( length( Data.( [planModeWord  '_Time' ] ) ), Config.Burst_NCells);
    K = 4;
end

Name = ['X','Y','Z'];
%ENU = zeros( K, Config.Burst_NCells);
xyz = zeros( K, Config.Burst_NCells);

for sampleIndex = 1:length(Data.( [planModeWord  '_Status' ]))
   if (bitand(bitshift(uint32(Data.( [planModeWord  '_Status' ])(sampleIndex)), -25),7) == 5)
      signXYZ=[1 -1 -1 -1];
   else
      signXYZ=[1 1 1 1];
   end

   hh = pi*(Data.([planModeWord  '_Heading_Cal'])(sampleIndex)+declination_correction-90)/180;
   pp = pi*Data.([planModeWord  '_Pitch'])(sampleIndex)/180;
   rr = pi*Data.([planModeWord  '_Roll'])(sampleIndex)/180;

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

   for ii = 1:K
      if (twoZs == 1) && (ii >= 3)
         axs = [ Name(3) num2str((ii-2),1) ];
      else
         axs = Name(ii);
      end
      xyz( ii, : ) = signXYZ(ii) * Data.( [ dataModeWord '_Vel' axs] )( sampleIndex, : )';
   end
   ENU = xyz2enu * xyz;
   EAllCells( sampleIndex, : ) = ENU( 1, : )';	
   NAllCells( sampleIndex, : ) = ENU( 2, : )';
   UAllCells( sampleIndex, : ) = ENU( 3, : )';
   if twoZs == 1
      U2AllCells( sampleIndex, : ) = ENU( 4, : )';
   end
end
%Config.( [ad2cpstr   dataModeWord '_CoordSystem' ] ) = 'enu';
Data.( [ dataModeWord '_VelEast' ] ) = EAllCells;
Data.( [ dataModeWord '_VelNorth' ] ) = NAllCells;
if twoZs == 1
	Data.( [ dataModeWord '_VelUp1' ] ) = UAllCells;
	Data.( [ dataModeWord '_VelUp2' ] ) = U2AllCells;
else
	Data.( [ dataModeWord '_VelUp' ] ) = UAllCells;
end
function [Data] = CalibrateCompass8_2017(Data, Par, planModeWord)
cor = 0;

% DataNew = Data;

HxHyHz = [294 -55 0]; % correction provided by Nortek based on past data
%from Jetyak. 
disp(['Nortek provided correction: HX=', num2str(HxHyHz(1)), ', HY=', num2str(HxHyHz(2)), ', HZ=', num2str(HxHyHz(3))])

% correction based on current data set
HxHyHz = [-Par(1) -Par(2) 0];
disp(['Data set specific correction: HX=', num2str(HxHyHz(1)), ', HY=', num2str(HxHyHz(2)), ', HZ=', num2str(HxHyHz(3))])

% disp(['SETUSER,HX=' num2str(HxHyHz(1)) ',HY=' num2str(HxHyHz(2)) ',HZ=' num2str(HxHyHz(3))])
% disp('SAVE,USER')

%HxHyHz = [301.5945 -77.7107 0]; % correction provided by Nortek based on past data

Data.( [ planModeWord '_Magnetometer_Cal' ] )(:,1) = Data.( [ planModeWord '_Magnetometer' ] )(:,1) + HxHyHz(1);
Data.( [ planModeWord '_Magnetometer_Cal' ] )(:,2) = Data.( [ planModeWord '_Magnetometer' ] )(:,2) + HxHyHz(2);
Data.( [ planModeWord '_Heading_Cal' ] ) = 180*atan2(Data.( [ planModeWord '_Magnetometer_Cal' ] )(:,2),Data.( [ planModeWord '_Magnetometer_Cal' ] )(:,1))/pi  + cor;
ii = find(Data.( [ planModeWord '_Heading_Cal' ] ) < 0);
Data.( [ planModeWord '_Heading_Cal' ] )(ii) = Data.( [ planModeWord '_Heading_Cal' ] )(ii) + 360;

% DataNew.BurstBT_MagnetometerX = DataNew.BurstBT_MagnetometerX + HxHyHz(1);
% DataNew.BurstBT_MagnetometerY = DataNew.BurstBT_MagnetometerY + HxHyHz(2);
% DataNew.BurstBT_Heading = 180*atan2(DataNew.BurstBT_MagnetometerY,DataNew.BurstBT_MagnetometerX)/pi + cor;
% i = find(DataNew.BurstBT_Heading < 0);
% DataNew.BurstBT_Heading(i) = DataNew.BurstBT_Heading(i) + 360;
% 
% DataNew.IBurst_MagnetometerX = DataNew.IBurst_MagnetometerX + HxHyHz(1);
% DataNew.IBurst_MagnetometerY = DataNew.IBurst_MagnetometerY + HxHyHz(2);
% DataNew.IBurst_Heading = 180*atan2(DataNew.IBurst_MagnetometerY,DataNew.IBurst_MagnetometerX)/pi + cor;
% i = find(DataNew.IBurst_Heading < 0);
% DataNew.IBurst_Heading(i) = DataNew.IBurst_Heading(i) + 360;

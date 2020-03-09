% name:         ADCP_postproc_PART2
% function:     2nd overall loop script
                %can be applied after indexing
                %misses only rot, but for the rest ready to make figures

close all
clear
clc

drpbx     = 'c:\users\wmkra\Dropbox\';
addpath([drpbx,'north_river\north_river_2017\Sig1000data\Functions']);
addpath([drpbx,'mrocky']); 
addpath([drpbx,'mwouter']);


for adcpday = 206;% [108, 109, 117, 137,138,205,206,209,212];%[108, 109]; %, %108
    disp(['adcpday=',num2str(adcpday)]);
    yrday = adcpday;
    daynr = yrday;
    
    PROC21_LoadStuff
    
%     PROC05_Arrows                     %figs only
%     PROC06_Rotation                   %proc (theta&rot)
%     PROC06_Rotation_Figs              %figs only
%     PROC07_VelocityFigs               %figs only
%     PROC07b_VelocityFigs_rotfixed     %figs only


%       PROC26_FigsGathered               %figs only, just new ones for the improved cases
      PROC27b_VelocityFigs_rotfixed     %figs only      
end

    
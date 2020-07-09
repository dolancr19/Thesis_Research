%% Prepare workspace
clear variables
clc

%% Initialize directories
%directory=["D:\Documents\Thesis_Research\Bigelow_10_30\sandshark_5_30_10_2019_____16_12_04_alvtmp\","D:\Documents\Thesis_Research\Bigelow_10_30\sandshark_5_30_10_2019_____16_22_41_alvtmp\","D:\Documents\Thesis_Research\Bigelow_10_30\sandshark_5_30_10_2019_____16_41_47_alvtmp\"];
%directory1='D:\Documents\Thesis_Research\Bigelow_10_30\';
% directory=["D:\Documents\Thesis_Research\Bigelow_11_26\sandshark_5_26_11_2019_____19_59_29_alvtmp\"];
% directory1='D:\Documents\Thesis_Research\Bigelow_11_26\';
% 12/12 Platypus
% directory=["D:\Documents\Thesis_Research\Bigelow_12_12\raw_data\PLATYPUS\sandshark_5_12_12_2019_____18_42_44_alvtmp\", "D:\Documents\Thesis_Research\Bigelow_12_12\raw_data\PLATYPUS\sandshark_5_12_12_2019_____19_36_40_alvtmp\", "D:\Documents\Thesis_Research\Bigelow_12_12\raw_data\PLATYPUS\sandshark_5_12_12_2019_____19_54_17_alvtmp\"];
% directory1='D:\Documents\Thesis_Research\Bigelow_12_12\processed_data\PLATYPUS\';
% % 12/12 Quokka
% directory=["D:\Documents\Thesis_Research\Bigelow_12_12\raw_data\QUOKKA\sandshark_6_12_12_2019_____18_51_04_alvtmp\", "D:\Documents\Thesis_Research\Bigelow_12_12\raw_data\QUOKKA\sandshark_6_12_12_2019_____19_09_45_alvtmp\", "D:\Documents\Thesis_Research\Bigelow_12_12\raw_data\QUOKKA\sandshark_6_12_12_2019_____19_16_39_alvtmp\"];
% directory1='D:\Documents\Thesis_Research\Bigelow_12_12\processed_data\QUOKKA\';
% 2/20 Platypus
% directory=["D:\Documents\Thesis_Research\Bigelow_2_20\raw_data\platypus\sandshark_5_20_2_2020_____16_39_21_alvtmp\","D:\Documents\Thesis_Research\Bigelow_2_20\raw_data\platypus\sandshark_5_20_2_2020_____16_43_05_alvtmp\","D:\Documents\Thesis_Research\Bigelow_2_20\raw_data\platypus\sandshark_5_20_2_2020_____16_48_48_alvtmp\","D:\Documents\Thesis_Research\Bigelow_2_20\raw_data\platypus\sandshark_5_20_2_2020_____17_08_50_alvtmp\"];
% directory1='D:\Documents\Thesis_Research\Bigelow_2_20\processed_data\platypus\';
% 2/20 Quokka
% directory=["D:\Documents\Thesis_Research\Bigelow_2_20\raw_data\quokka\sandshark_6_20_2_2020_____17_01_58\"];
% directory1='D:\Documents\Thesis_Research\Bigelow_2_20\processed_data\quokka\';
% MIT Platypus
% directory=["D:\Documents\Thesis_Research\MIT_Sailing\raw_data\platypus\sandshark_5_14_9_2018_____13_47_51_unix_time_unix_time_gps_spliced_alvtmp\"];
% directory1='D:\Documents\Thesis_Research\MIT_Sailing\processed_data\platypus\';
% MIT Quokka
% directory=["D:\Documents\Thesis_Research\MIT_Sailing\raw_data\quokka\sandshark_6_21_9_2018_____14_05_08_unix_time_alvtmp\"];
% directory1='D:\Documents\Thesis_Research\MIT_Sailing\processed_data\quokka\';
% directory=["D:\Documents\Thesis_Research\MIT_Sailing\raw_data\quokka\sandshark_6_14_9_2018_____13_36_14_unix_time_alvtmp\"];
% directory1='D:\Documents\Thesis_Research\MIT_Sailing\processed_data\quokka\';
directory=["D:\Documents\Thesis_Research\MIT_Sailing\raw_data\wombat\sandshark_7_14_9_2018_____13_49_18_unix_time_alvtmp\"];
directory1='D:\Documents\Thesis_Research\MIT_Sailing\processed_data\wombat\';

%% Specify files to import
moos_variable=["NAV_SPEED","NAV_HEADING","NAV_X","NAV_Y","NAV_LAT","NAV_LONG","ACOUSTIC_BEARING","ACOUSTIC_RANGE","ACOUSTIC_SOURCE_X","ACOUSTIC_SOURCE_Y", "NAV_DEPTH","ACOUSTIC_STDDEV_1","ACOUSTIC_STDDEV_2"];
%start_time=[1572451923.49,1572452560.78,1572453706.79];
%start_time=[1574798368.75]; %11-26
%start_time=[1576176163.79,1576179399.92,1576180456.34]; %12-12 Platypus
%start_time=[1576176663.37, 1576177784.93, 1576178198.67]; %12-12 Quokka
% start_time=[1582216760.88,1582216984.61,1582217328.16,1582218529.62]; %2-20 Platypus
% start_time=[1582218117.55]; %2-20 quokka
% start_time=[1536932870.64]; %MIT Platypus
% start_time=[1537538708.17]; %MIT Quokka
% start_time=[1536932173.36]; %MIT Quokka 9/14/18
start_time=[1536932958.13]; %MIT Wombat 9/14/18


% Right now, the following format is hard coded for the type of files we are
    % using.  It can be made more flexible by adding variables to describe
    % the column names and data types.
    
%% Generate .mat files for each specified variable
for jj=1:length(directory)
    for ii=1:length(moos_variable)
        filename=directory(jj) + moos_variable(ii) + '.klog';

        var = convertStringsToChars(moos_variable(ii));
            varNames = {'TIME','variable','source',var};
            varTypes = {'double','char','char','double'};
            delimiter = ' ';

        dataStartLine = 1;
        extraColRule = 'ignore';
        ConsecutiveDelimitersRule='join';

        opts = delimitedTextImportOptions('VariableNames',varNames,...
                                        'VariableTypes',varTypes,...
                                        'Delimiter',delimiter,...
                                        'DataLines', dataStartLine,...
                                        'ExtraColumnsRule',extraColRule, 'ConsecutiveDelimitersRule',ConsecutiveDelimitersRule);
        working=readtable(filename,opts);
        if ii==1 || ii==3 || ii==5 || ii==7 || ii==9 || ii==12
            data=removevars(working,{'variable','source'});
        elseif ii==2
            working=removevars(working,{'TIME','variable','source'});
            data=[data working];
%             data.TIME=data.TIME+start_time(jj);
            save(directory1 + string(start_time(jj)) + '_NAVDR.mat','data','-v7.3')
        elseif ii==4
            working=removevars(working,{'TIME','variable','source'});
            data=[data working];
%             data.TIME=data.TIME+start_time(jj);
            save(directory1 + string(start_time(jj)) + '_NAVXY.mat','data','-v7.3')
        elseif ii==6
            working=removevars(working,{'TIME','variable','source'});
            data=[data working];
%             data.TIME=data.TIME+start_time(jj);
            save(directory1 + string(start_time(jj)) + '_NAVLL.mat','data','-v7.3')
            s=geoshape(data.NAV_LAT,data.NAV_LONG);
            filename_kml=directory1 + string(start_time(jj)) + "_NAV.kml";
            kmlwrite(filename_kml,s);
        elseif ii==8
            working=removevars(working,{'TIME','variable','source'});
            data=[data working];
%             data.TIME=data.TIME+start_time(jj);
            save(directory1 + string(start_time(jj)) + '_ACOUSTICRB.mat','data','-v7.3')
        elseif ii==10
            working=removevars(working,{'TIME','variable','source'});
            data=[data working];
%             data.TIME=data.TIME+start_time(jj);
            save(directory1 + string(start_time(jj)) + '_SOURCEXY.mat','data','-v7.3')
        elseif ii==11
            data=removevars(working,{'variable','source'});
%             data.TIME=data.TIME+start_time(jj);
            save(directory1 + string(start_time(jj)) + '_NAVDEPTH.mat','data','-v7.3')
        elseif ii==13
            working=removevars(working,{'TIME','variable','source'});
%             working(1,:)=[]; %only if row numbers not equal
            data=[data working];
%             data.TIME=data.TIME+start_time(jj);
            save(directory1 + string(start_time(jj)) + '_ACOUSTIC_STDDEV.mat','data','-v7.3')
        end
    end
end
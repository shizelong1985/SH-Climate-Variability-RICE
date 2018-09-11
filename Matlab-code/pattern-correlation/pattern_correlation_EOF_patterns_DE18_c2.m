%%%%%%%%%%%%%%
% pattern_correlation_EOF_patterns_DE18_c2.m
% Daniel Emanuelsson
% Matlab 2017a
% Github version 1  11-09-2018
%%%%%%%%%%%%%%
clear all
close all
%%%%
%% Save patterns from ERA and HadISST regression
% and HadISST regression files

Param_nr=1; % (1) SAT (2) SIC

if Param_nr==1 % SAT
        savefilename = ['C:\PHD\matlab_storage_of_output_files\corr_patterns_',filename,'_c.mat']; % file name from regression file, run this file first
        save(savefilename,'s_all_c','era_lat','era_long');
        
elseif Param_nr==2 % SIC      
        savefilename = ['C:\PHD\matlab_storage_of_output_files\corr_patterns_',filename,'_c.mat']; % file name from regression file, run this file first
        save(savefilename,'s_all','HadISST_lat','HadISST_lon');   

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load saved patterns
% compare pattern pairs in Fig. 3 and Fig. 7
Param_nr=1; % (1) SAT (2) SIC
PC_nr=3; % (1), (2), (3) %%%%%%%%%%%%%


if Param_nr==1 % SAT
%%%%%%%%%%%%%%%%%%%%%%%%%    
% for SAT
%%%%%%%%%%%%%%%%%%%%%%%%%
% load saved surfaces
% Note the 1, 3, 2 order. Due to the SAT patterns relation to z500
% surface #1
    if PC_nr==1 % SAT SAM regression surface
        load('C:\PHD\matlab_storage_of_output_files\corr_patterns_ERA_interim_regression_2mT_SAM_Annual_Monthly_c.mat')
    elseif PC_nr==3 % SAT PSA1
        load('C:\PHD\matlab_storage_of_output_files\corr_patterns_ERA_interim_regression_2mT_PSA1_Annual_Monthly_c.mat')
    elseif PC_nr==2 % SAT PSA2
        load('C:\PHD\matlab_storage_of_output_files\corr_patterns_ERA_interim_regression_2mT_PSA2_Annual_Monthly_c.mat')
    end

    A_s_all_c=s_all_c;

    lat_nr=2; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if lat_nr==1
        A_s_all_c=A_s_all_c(:,(161:end)); % 30S
    elseif lat_nr==2
        A_s_all_c=A_s_all_c(:,(181:end)); % 45S
    end

    %%%%%%%%%%%%%%%%%
    % surface #2
    if PC_nr==1 % SAT EOF1 surface
        load('C:\PHD\matlab_storage_of_output_files\corr_patterns_ERA_interim_regression_2mT_PC1_annual_Monthly_c.mat')
    elseif PC_nr==2 % SAT EOF2
        load('C:\PHD\matlab_storage_of_output_files\corr_patterns_ERA_interim_regression_2mT_PC2_annual_Monthly_c.mat')
    elseif PC_nr==3 % SAT EOF3
        load('C:\PHD\matlab_storage_of_output_files\corr_patterns_ERA_interim_regression_2mT_PC3_annual_Monthly_c.mat')
    end

    B_s_all_c=s_all_c;


    if lat_nr==1
        B_s_all_c=B_s_all_c(:,(161:end)); % south of 30S
    elseif lat_nr==2
        B_s_all_c=B_s_all_c(:,(181:end)); % south of 45S    
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for SIC
elseif Param_nr==2  % SIC
    
    %%%%%%%%%%%%%%%%%%
    % surface #1
    if PC_nr==1 % SIC SAM regression surface
        load('C:\PHD\matlab_storage_of_output_files\corr_patterns_HadISST_regression_SIC_SAM__monthly_1979-2014_c.mat');
    elseif PC_nr==2 % SIC PSA1  
        load('C:\PHD\matlab_storage_of_output_files\corr_patterns_HadISST_regression_SIC_PSA1__monthly_1979-2014_c.mat'); 
     elseif PC_nr==3 % SIC PSA2    
         load('C:\PHD\matlab_storage_of_output_files\corr_patterns_HadISST_regression_SIC_PSA2__monthly_1979-2014_c.mat');
    end
    
    A_s_all_c=s_all;

    A_s_all_c=A_s_all_c(:,(151:end)); % south of 60S
    
    
    %%%%%%%%%%%%%%%%%%
    % surface #2   
    if PC_nr==1      % SIC EOF1 surface
        load('C:\PHD\matlab_storage_of_output_files\corr_patterns_HadISST_regression_SIC_PC1__monthly_1979-2014_c.mat');
    elseif PC_nr==2    % SIC EOF2 
        load('C:\PHD\matlab_storage_of_output_files\corr_patterns_HadISST_regression_SIC_PC2__monthly_1979-2014_c.mat');
     elseif PC_nr==3    % SIC EOF3
        load('C:\PHD\matlab_storage_of_output_files\corr_patterns_HadISST_regression_SIC_PC3__monthly_1979-2014_c.mat');        
    end
    
    B_s_all_c=s_all;
    B_s_all_c=B_s_all_c(:,(151:end)); % south of 60S

end

%%%%%%%%%%%%%%%%%%%%%%%%
% Y=xcorr2(A_s_all_c,B_s_all_c);
% Y=corr2(A_s_all_c, B_s_all_c);

Y=corr2_cDE(A_s_all_c, B_s_all_c)     
% replaced mean2(a) with nanmean(nanmean(a)) otherwise everthing else is
% the same as corr2

% edit corr2
% edit corr2_cDE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
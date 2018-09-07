%%%%%%%%%%%%%%%%%
%  eof_code_run_annual_seasonal_DE18_c1.m
%  EOFs SH
%  Daniel Emanuelsson
%  Matlab 2017a
%  Github version
%%%%%%%%%
% subfunctions
% * keep_var.m        UoW online archive Atmospheric Science
% * annave.m          UoW online archive Atmospheric Science
% * eof_j.m           James R
% * annual_mean_DE.m  DE
% * HadISST_load_lat_lon.m DE
%%%%%%%%%%%%%%%%%
clear all
close all

%%%%%%%%%%%%%%%%%%%%
sea_nr=1;   % all values are used here
%%%%%%%%%%%%%%%%%%%%%%
reanalysis_nr=4; % (1) * ERA-interim (3) * HadISST SIC

    if reanalysis_nr==1
        reanalysis_str='ERA-Interim';
    elseif reanalysis_nr==3
        reanalysis_str='HadISST';
    end
yr_s=1979;
yr_e=2014;%%%%%%%%%%%%%%%
        
%     if sea_nr==1
    season='annual'; 
%     end
    
%%%%% ERA-Interim %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 addpath C:\PHD\ERA_interim
 
if reanalysis_nr==1 
    era_time=ncread('ERA_int_monthly_z500_2.nc','time');
    era_long=ncread('ERA_int_monthly_z500_2.nc','longitude');
    era_lat=ncread('ERA_int_monthly_z500_2.nc','latitude');
 
    dayssince111=era_time/24;
    datevalue=dayssince111+datenum(1900,1,1);
    date_vec=datevec(double(datevalue)); 
    yyyy=date_vec(:,1); 
    mm=date_vec(:,2);
    era_year_num=yyyy+(mm/12-1/12);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif reanalysis_nr==3 % HadISST
    % edit HadISST_load_lat_lon
    % nameoffile='HadISST_ice_c.nc';
    nameoffile='HadISST_ice_c2.nc';
    [HadISST_lon, HadISST_lat, HadISST_time, HadISST_year_num, mm]=HadISST_load_lat_lon(nameoffile);
   % use excisting names
    era_long=HadISST_lon;
    era_lat=HadISST_lat;
    era_time=HadISST_time; 
    era_year_num=HadISST_year_num;  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

era_count=find(era_year_num==yr_e)+11;  % data in months
era_start=find(era_year_num==yr_s); 

mm_in=mm(era_start:era_count);  % seasonal index
mm_in_c=[era_start:era_count]';
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
geo_lev=1; % (1) z500 (7) 2mT (8) SIC 

if geo_lev==1
    geo_str='z500';
elseif geo_lev==7 % 2mT SAT  
    geo_str='2mT';
elseif geo_lev==8 % SIC  
    geo_str='SIC';       
end

name=geo_str;
 
if  geo_lev<=4 && reanalysis_nr==1    
        era_z500=ncread(strcat('ERA_int_monthly_',geo_str,'_2.nc'),'z'); 
      
elseif geo_lev==7
        era_z500=ncread(strcat('ERA_int_monthly_',geo_str,'_2.nc'),'t2m'); % added 2017
           
elseif geo_lev==8 % SIC
        era_z500=ncread('HadISST_ice_c2.nc','sic'); % added 2017
        
        X=era_z500;
        
        [A B C]=size(X);
            for i=1:A
                for j=1:B
        %%%%%%%%%%%%%%%%%%
        %%% Replace NaNs
                    if  isnan(X(i,j,:))==1
                        X(i,j,:)=0;
                    end
        %%%%%%%%%%%%%%%%%%%%
                end
            end
        
    era_z500=X;    
        
end
 
    if geo_lev==1 &&  reanalysis_nr==1    
        era_z500=era_z500/9.80665;
    end
   
 era_z500=era_z500(:,:,(era_start:era_count));
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%
%     [xkeep, ykeep] = keep_var(lim, x, y);
%    where lim = [minx maxx miny maxy] to be kept.
if geo_lev==1
  lims = [0 360 -90 -20]; % z500

elseif geo_lev==7
  lims = [0 360 -90 -30]; % Yu et al 2015 2mt

elseif geo_lev==8
   lims = [-150 -30 -75 -64]; % sic---in use
end

[xk, yk] = keep_var(lims, era_long, era_lat);
era_z500 = era_z500(xk, yk, :);
era_lat = era_lat(yk); 
era_long =era_long(xk);

% era_z500(:,(1:147),:)=[];  % Only use the SH <20S
% era_lat(1:147)=[]; 
% 
% %%%%%
% era_z500(:,(69:94),:)=[];  % Remove 70 to 90S   %%%% Boundaries from Marshall 2016
% era_lat(69:94)=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% resize grid
 
 if reanalysis_nr==1  % Only for ERA-Interim (high resolution data)
 
    [a b c]=size(era_z500);
 
    scale=0.5;
    %scale=1;
 
        for i=1:c
        
            era_z500_v(:,:,i)=resizem(era_z500(:,:,i),scale);
    
        end 
 
    [av, bv, cv]=size(era_z500_v);
    xq = (1:bv)*1/scale; % Lat
    yq = (1:av)*1/scale; % Lon

    % resize lat long
    x=[1:length(era_lat)]';  % Lat
    [p,s]=polyfit(x,era_lat,1);

    era_lat_v=polyval(p,xq');
    clear x p s

    x=[1:length(era_long)]';
    [p,s]=polyfit(x,era_long,1);

    era_lon_v=polyval(p,yq');    
%     clear a b
    
    era_z500=era_z500_v;
    era_lat=era_lat_v;
    era_long=era_lon_v;  
 end
 
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
era_lat_c=era_lat;
[nlon, nlat, ntim]= size(era_z500);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Order of these three steps doesn't matter
%%% -=(1)=-
%  z = cosweight(z, lat);
% weighted to account for the decrease in area towards the pole.

%  era_z500_c = reshape(era_z500, ntim, nlat, nlon);  %  z500_eof

era_z500_c = permute(era_z500, [3 2 1]);  %  z500_eof

%%%%%%%%%%%%%%%
%     B = permute(A,ORDER) rearranges the dimensions of A so that they
%     are in the order specified by the vector ORDER.  The array produced
%     has the same values as A but the order of the subscripts needed to 
%     access any particular element are rearranged as specified by ORDER.
%%%%%%%%%%%%%%%%%

era_z500_c2=cosweight(era_z500_c,era_lat_c); % Original UoW function    (time x lat x lon) 

%%%%%%%%%%%era_z500_c2=cosweight_c(era_z500,era_lat_c);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
 %%% -=(2)=-
% Remove monthly climatology for all cells 

 era_z500_c3 = double(reshape(era_z500_c2, ntim, nlat*nlon));  % One time series (column) for each grid point
 
 %  lat 1, 2, 3, 4,..........33...lat 1, 2, 3, 4.....33       until 33X240=7920
 %  
 %1
 %time
 %
 %444
 
 [era_z500_c4,clim_z500] = annave(era_z500_c3);   % checked, Removes
%  remove seasonal cycle

%%%%%%%%%% SAVE used in regression -= do not downscale above =-%%%%%%%%%%%%%%%%%% 
%  filedir ='C:\PHD\matlab_storage_of_output_files\';
%  name_c=['clim_',name,'_',num2str(yr_s),'_',num2str(yr_e)];
%  name_c=['clim_',name,'_',reanalysis_str,'_',num2str(yr_s),'_',num2str(yr_e)];

% savefilename =[filedir, name_c, '.mat'];
% save(savefilename,'era_z500_c4','clim_z500'); 

% stop here for regression save
%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%  Don't detrend if you want to keep trend.....
 
 
 detrend_nr=0; % (1) detrend (0) keep trend
 
 if  detrend_nr==0
    if sea_nr==1 
    era_z500_c8=   era_z500_c4;
    elseif sea_nr>=2 
    era_z500_c5=   era_z500_c4;  
    end
    
 elseif detrend_nr==1
 %%% -=(3)=-                   *    Order of 2 and 3 doesn't seem to matter
    if sea_nr==1
    era_z500_c8= detrend(era_z500_c4); % Linear detrend
    else % seasonal
    era_z500_c5= detrend(era_z500_c4); % Linear detrend   
    end
 end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%%
%%%%%%  James EOF function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 eof_out=eof_j(era_z500_c8, 10);% ,[],'svd'); 

 % edit eof_j
% out = eof(ldata,nsave,scaling,atype)
%
% Get eigenvectors/values by eig or SVD method
% Input:
%	ldata	-	n by p data array
%	nsave	-	How many vectors (etc) to save/return - 0 => all
%	scaling -	if true (1), scale columns of ldata by std(ldata) first
%	atype	-	analysis type: 'eig' or 'svd' (default depends on n/p)
%
% Output:
%  Structure with components
%	vect		-	Eigenvectors (Most important "nsave")
%	val		-	Eigenvalues (  ..     ..       ..   )
%	series	-	Principal components (loadings, time series)
%	regmap	-	Regression map version of vectors
%	totvar	-	Total variance of input data
%	expv		-	Percent explained variance
%	rank		-	Rank info: [sum(val>1e-7) spatial d.o.f]
%	std		-	standard deviations of columns of ldata (=1 if scaling)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%
  
MA_PC_time_series=eof_out.series;
era_year_num=era_year_num(mm_in_c); %%% start from the beginning if you need to run it again


%%%%%%%%%% standardize using zscore
MA_PC_time_series_standardized=zscore(MA_PC_time_series);  
% MA_PC_time_series_standardized=MA_PC_time_series_standardized(mm_in_c,:);

% edit annual_mean_DE

MA_PCs_save=annual_mean_DE(yr_s, yr_e,era_year_num, MA_PC_time_series_standardized,sea_nr );


clear a b c
[a, b, c]=size(MA_PC_time_series_standardized);

MA_PCs_save_monthly=zeros(a,b+1);
MA_PCs_save_monthly(:,1)=era_year_num;
MA_PCs_save_monthly(:,(2:11))=MA_PC_time_series_standardized; % save for 10 fist EOFs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Save Data 

 filedir ='C:\PHD\matlab_storage_of_output_files\';
 name_end_str='_c2';
 
 name_str=[reanalysis_str,'_PCs','_', geo_str,'_lim',num2str(lims(1)),...
     '-',num2str(lims(2)),'_',num2str(lims(4)),'_',num2str(lims(3)),'_',num2str(yr_s),'-',num2str(yr_e),'_',season,name_end_str];

    if  detrend_nr==1
        savefilename =[filedir, name_str '.mat'];
    elseif  detrend_nr==0
        savefilename =[filedir, name_str,'_trend' ,'.mat'];     
    end
    
                save(savefilename,'MA_PCs_save','MA_PCs_save_monthly');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%
%  HadISST_load_lat_lon.m
%  Daniel Emanuelsson
%  Matlab 2017a
%  Github version 1
%  Loads lat lon time HadISST v. 1.1 data 
%%%%%%%%%%%%%%%%%%%%%%%

function  [HadISST_lon, HadISST_lat, HadISST_time, HadISST_year_num, mm]=HadISST_load_lat_lon(nameoffile)

addpath C:\PHD\HadISST_data\ % Ienovo
addpath C:\PHD\matlab_lib\Data\

%   ncdisp('HadISST_ice_c.nc');
   
 HadISST_time=ncread('HadISST_ice_c.nc','time'); %units         = 'days since 1870-1-1 0:0:0'
 % missing values  missing_value = -1e+30 ( looks more like -1000)\
 %monthly values
 
 HadISST_time_c=HadISST_time+datenum(1870,1,1);
 
 date_vec=datevec(double(HadISST_time_c)); 
 yyyy=date_vec(:,1); 
 mm=date_vec(:,2);
 HadISST_year_num=yyyy+(mm/12-1/12); 

 HadISST_time_bnds =ncread('HadISST_ice_c.nc','time_bnds');
 HadISST_lat=ncread('HadISST_ice_c.nc','latitude');
 HadISST_lon=ncread('HadISST_ice_c.nc','longitude');
end
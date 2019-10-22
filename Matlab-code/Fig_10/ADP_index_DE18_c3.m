%%%%%%%%%%%%%%%%%%%%%%
% ADP_index_DE18_c3.m
% Daniel Emanuelsson
% Matlab 2018a
% Github version 1 [uploaded]
%%%%%%%%%%%%%%%%%%%%%
% Antarctic Dipole (ADP)
% Ref paper: Yuan 2004
% ENSO-related impacts on Antarctic sea ice: a synthesis of
% phenomenon and mechanisms
% Antarctic Science 16 (4): 415–425 (2004)
%%%%%%%%%%%%%%
% subfunctions
%
%* keep_var.m            UoW online archive Atmospheric Science
%* annual_mean_DE.m      DE
%* p_level.m             DE
%* corrcoef_df.m   UoW Steig
%%%%%%%%%%%%%%%
clear all
close all

param_nr=2;
if param_nr==1
    name='SST';
elseif param_nr==2
    name='SIC';
end

%%%%%%%%%%%%%%%%%
 yr_s=1979;
 yr_e=2014;
 %%%%%%%%%%%%%%%%%%%
 
addpath C:\PHD\HadISST_data\ % Ienovo

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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if strcmp(name,'SST')==1
            %%%%%%%%%% file    _c contains data until 2015
            %name='SST';
            HadISST_sst=ncread('HadISST_sst_c.nc','sst'); % 360x180x1728 lon x lat x time
            M=HadISST_sst;
        elseif strcmp(name,'SIC')==1
            %   name='ice'; 
            HadISST_ice=ncread('HadISST_ice_c.nc','sic'); % 360x180x1728 lon x lat x time
            M=HadISST_ice;
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%%
%     [xkeep, ykeep] = keep_var(lim, x, y);
%    where lim = [minx maxx miny maxy] to be kept.
 
lims = [-150 -130 -70 -60]; % 1.  SIC ADP Pacific center 

[xk, yk] = keep_var(lims, HadISST_lon, HadISST_lat);
M_c =M(xk, yk, :);
HadISST_lat_c = HadISST_lat(yk); 
HadISST_lon_c =HadISST_lon(xk); 
       
SIC_mean_Pacific=squeeze(nanmean(nanmean(M_c)));
%%%%%%%%%%
lims2 = [-60 -40 -70 -60]; % 2.  SIC ADP Atlantic center 
  
[xk, yk] = keep_var(lims2, HadISST_lon, HadISST_lat);
M_c2 =M(xk, yk, :);
HadISST_lat_c2 = HadISST_lat(yk); 
HadISST_lon_c2 =HadISST_lon(xk); 

SIC_mean_Atlantic=squeeze(nanmean(nanmean(M_c2)));

ADP_index=SIC_mean_Pacific-SIC_mean_Atlantic;

% HadISST_year_num
HadISST_count=find(HadISST_year_num==yr_e)+11; % values in months
% HadISST_start=289; %1894
HadISST_start=find(HadISST_year_num==yr_s); % 

ADP_index_c=ADP_index(HadISST_start:HadISST_count);
HadISST_year_num_c=HadISST_year_num(HadISST_start:HadISST_count);

% help annual_mean_DE.m
% edit annual_mean_DE.m    
sea_nr=1; % annual
MA_ADP_save=annual_mean_DE(yr_s, yr_e,HadISST_year_num, ADP_index,sea_nr);
%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure
% check data
  h=fig('units','inches','width',14,'height',5,'font','Helvetica','fontsize',16,'border','on');
     h2=plot(MA_ADP_save(:,1), MA_ADP_save(:,2),'.-','LineWidth',4,'MarkerSize',14,'Color',[0.9,0.2,0.0]); % annual
     hold on
     h3=plot(HadISST_year_num, ADP_index,'.-','LineWidth',2,'MarkerSize',14,'Color',[0.2,0.2,0.2]); % monthly
     set(gca,'XLim',[1977 2012])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ADP Moving Variance Fig.

 moving_window=11; % years
 tim_cp_nr=1; % (1) Annual (2) Monthly
 
 if tim_cp_nr==2
    v_cp_c= movingvar(ADP_index_c,moving_window*12);  
 elseif tim_cp_nr==1    
%   v_cp_c= movingvar( MA_ADP_save(:,2),moving_window); % annual
    v_cp_c= movvar(MA_ADP_save(:,2),moving_window); % annual
 end

 h=fig('units','inches','width',14,'height',5,'font','Helvetica','fontsize',16,'border','on');
 hold on
 
            vl= vline([1994, 2004]);
            set(vl,'LineWidth',3,'color',[0.3 0.5 1],'LineStyle','--');
 
  if tim_cp_nr==2
        h2=plot(HadISST_year_num_c, v_cp_c,'.-','LineWidth',4,'MarkerSize',14,'Color',[0.9,0.2,0.0]);      
  elseif tim_cp_nr==1
%        h2=plot( MA_ADP_save(:,1), v_cp_c,'.-','LineWidth',4,'MarkerSize',14,'Color',[0.9,0.2,0.0]);
         h2=plot( MA_ADP_save((6:30),1), v_cp_c(6:30),'.-','LineWidth',4,'MarkerSize',14,'Color',[0.9,0.2,0.0]); % cut out ends with the same values
  end
 
         set(gca,...
          'linewidth',3,...
          'FontWeight','bold' );
 
%        xlabel('Year','FontWeight','bold','FontSize',18 );
         ylabel('Variance','FontWeight','bold','FontSize',20 );
         box on
 %%%%%%%%%%%%%%%%%%%%%%%
 % letter      
  letter_pos=[65 380 40 60];
   letter='';
          TextBox = uicontrol('style','text');
          set(TextBox,'String',letter,'position',letter_pos,'FontWeight','bold','FontSize',34 ); % x position y position xsize ysize
          set(TextBox,'foregroundcolor', [0 0 0], ...
         'backgroundcolor', [1 1 1]);
%%%%%%%%%%%%%%%%%%%%%%%%%    
% label 
 title_str='SIC ADP Var';
 axestext_c(0.9999999,0.012,title_str,'FontWeight','bold','FontSize',20 );
 %%%%%%%%%%%%%%%%%%%%%%%%%
           set(gca,'XLim',[1977 2012])
           set(gca,'YLim',[0.00 0.022])
%%%%%%%%%%%%%%%%%%%%%%%%%
% save fig
filename=['ADP_variance_w_',num2str(moving_window)];
filedir ='C:\PHD\matlab_storage_of_output_files\figures\';
savefilename_c=strcat(filedir,filename);
% save as png 
orient landscape
% increase quality of picture by increase to -r300
export_fig('-png','-nocrop','-painters', '-depsc','-opengl', '-r150',savefilename_c); % PNG-nocrop'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%% save index
save_ind=0;
if save_ind==1
        savefilename =strcat( 'C:\PHD\matlab_storage_of_output_files\Ma_ADP_index_HadISST_SIC',num2str(yr_s),'_',num2str(yr_e),'.mat'); 
        save(savefilename,'MA_ADP_save');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 

X=detrend(MA_ADP_save((1:33),2));

load('C:\PHD\matlab_storage_of_output_files\RICE_combined_Deep_1213B_c19.mat'); % Winstrup 2017 May age-scale
Y=detrend(MA_save((80:112),2));

%correlation
%  factor_nr
[r,p]=corrcoef_df(X,Y);
rp=[r(2),p(2)]
p_level(p(2))
%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%
% SAT_SIC_z500_RICE_DE18_c3.m
% dD corr w. SAT, SIC and z500
% Daniel Emanuelsson
% Matlab 2017a
% Github version 1 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions
% * IPO_vline.m DE
% * movingCorrelation_c.m DE
% * fig.m fileexchange
% * export_fig.m fileexchange
%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
%%%%%
addpath C:\PHD\matlab_lib\Data\      
 %%%%%%%%%%%%%% 
 % SAT ERA-Interim 
 %%%%%%%%%%%%%%%%%
    % ECMWF ERA-Inerim (Dee et al. 2011)
    % Monthly means of daily means
    % surface 2mT         
era_time=ncread('ERA_int_monthly_2mT_2.nc','time');
era_long=ncread('ERA_int_monthly_2mT_2.nc','longitude');
era_lat=ncread('ERA_int_monthly_2mT_2.nc','latitude');
 
dayssince111=era_time/24;
datevalue=dayssince111+datenum(1900,1,1);
date_vec=datevec(double(datevalue)); 
yyyy=date_vec(:,1); 
mm=date_vec(:,2);
era_year_num=yyyy+(mm/12-1/12);
geo_str='2mT';
      
SAT_ERA_Interim=ncread(strcat('ERA_int_monthly_',geo_str,'_2.nc'),'t2m'); % added 2017    

%%%%%%%%%%%RICE
SAT_ERA_Interim_c=squeeze(SAT_ERA_Interim(265,227,:));
% 456/12
SAT_ERA_Interim_c=SAT_ERA_Interim_c-273.15;

         SAT_ERA_Interim_c2=reshape(SAT_ERA_Interim_c,12,456/12);
         SAT_ERA_Interim_c3=permute(SAT_ERA_Interim_c2,[2 1]);
         SAT_ERA_Interim_annual= nanmean(SAT_ERA_Interim_c3,2);
 %%%%%%%% EA
         SAT_ERA_Interim_EA_c=squeeze(SAT_ERA_Interim(161,214,:));
        SAT_ERA_Interim_EA_c=SAT_ERA_Interim_EA_c-273.15;
         SAT_ERA_Interim_EA_c2=reshape(SAT_ERA_Interim_EA_c,12,456/12);
         SAT_ERA_Interim_EA_c3=permute(SAT_ERA_Interim_EA_c2,[2 1]);
         SAT_ERA_Interim_EA_annual= nanmean(SAT_ERA_Interim_EA_c3,2);
         
         ERA_Interim_time_annual=[1979:2016];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%% Figure SAT
% Fig 6b
 close all     
 h=fig('units','inches','width',14,'height',5,'font','Helvetica','fontsize',16,'border','on');  
   box on
   hold on      
            
       set(gca,'YLim',[-6 6])  
      wr_c=SAT_ERA_Interim_annual-nanmean(SAT_ERA_Interim_annual(1:31));
      h3=plot(ERA_Interim_time_annual, wr_c ,'-*','LineWidth',3,'MarkerSize',14,'Color',[0.0,0.2,0.7]); 
            
            addpath C:\PHD\matlab_lib\Library
            IPO_vline(10)
      
%%%%%%%%-=ISO=-%%%%%%%%%%%%%
add_yr=2; % (0) 2009 (2) 2011
%%%%Load RICE
RICE_nr=1;
 
    load('C:\PHD\matlab_storage_of_output_files\RICE_combined_Deep_1213B_c19.mat'); % May 2017
    stacked_record_annual_Ma=MA_save;
    iso_str='stack';
    
date_annual=stacked_record_annual_Ma((1:110+add_yr),1); 
X_in=stacked_record_annual_Ma((1:110+add_yr),2); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  RICE iso
 h16=plot(stacked_record_annual_Ma((1:110+add_yr),1),anomaly(X_in(1:end)),'*-','Color',[0.4,0.4,0.4],'LineWidth',3,'MarkerSize',14);
 set(gca,'XLim',[1977 2012])
%         set(gca,'YLim',[-30 -15])
        set(gca,...
          'linewidth',3,...
          'FontWeight','bold' );
      
%   ylabel('SAT Anom. R. 1979-2009','FontWeight','bold','FontSize',16 ); %  (°C)
  ylabel('Anomaly','FontWeight','bold','FontSize',20 ); %  (°C) 
  title_str='SAT';
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% label   
  letter_pos=[1125 447 100 40];
  
          TextBox = uicontrol('style','text');
          set(TextBox,'String',title_str,'position',letter_pos,'FontWeight','bold','FontSize',26 ); % x position y position xsize ysize
          set(TextBox,'foregroundcolor', [0 0 0], ...
         'backgroundcolor', [1 1 1]);     
     
%%%%%%%%%%%%%%%%%%%%%%%%%
hl=hline(0);
set(hl,'LineWidth',3); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % letter
 letter_pos=[50 380 60 60];  
 letter='b';
          TextBox = uicontrol('style','text');
          set(TextBox,'String',letter,'position',letter_pos,'FontWeight','bold','FontSize',34 ); % x position y position xsize ysize
          set(TextBox,'foregroundcolor', [0 0 0], ...
         'backgroundcolor', [1 1 1]);    
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % legend
     hs3='SAT ERA-I annual';
     hs4='SAT ERA-I 10-yr RM';
     hs16='{\delta}D annual';
     
h_leg=legend([  h3 h16 ],   hs3, hs16); 
set(h_leg,'Position',[0.27 0.80 0.1 0.1],'EdgeColor',[0 0 0],'FontSize',18 ); % ,'color','none'
    
 %%% Save Fig.
     filename=['SAT_RICE_anomaly','_', iso_str];  
     filedir ='C:\PHD\matlab_storage_of_output_files\figures\';
     savefilename_c=strcat(filedir,filename);

% save as png 
orient landscape
% increase quality of picture by increase to -r300
export_fig('-png','-nocrop','-painters', '-depsc','-opengl', '-r170',savefilename_c); % PNG-nocrop'  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fig 6c(c)   %%% moving correlation with RICE dD
sliding_w=11;%%%%%%%%%%%%%%%%%%%%%%%%%%%
add_yr=2; % (0) 2009 (2) 2011
    
load('C:\PHD\matlab_storage_of_output_files\ERA-Interim_PCs_2mT_lim0-360_-30_-90_1979-2014_annual.mat'); % PCs SAT ERA-Interim
MA_PCs_save_SAT=MA_PCs_save;
MA_PCs_save_monthly_SAT=MA_PCs_save_monthly;
clear MA_PCs_save MA_PCs_save_monthly

le=length(MA_PCs_save_SAT((1:33),1));

h=fig('units','inches','width',14,'height',5,'font','Helvetica','fontsize',16,'border','on');
   box on
   hold on   
   set(gca,'YLim',[-1 1])
 %%%%%%%%%%%%%%            
addpath C:\PHD\matlab_lib\Library
%   edit IPO_vline
   IPO_vline(10)
   hold on          
      %%%% ERA-Interim
      %%%%%%%%%%%%%%%%%%%%%%% RICE
le=length(stacked_record_annual_Ma(80:110+add_yr));

Ma_corr2=NaN(le,3); % use NaN instead of zero because NPGO index is shorter

% iso
Ma_corr2(:,1)=detrend(stacked_record_annual_Ma((80:110+add_yr),2));
Ma_corr2(:,2)=detrend(SAT_ERA_Interim_annual(1:31+add_yr));
Ma_corr2(:,3)=detrend(SAT_ERA_Interim_EA_annual(1:31+add_yr));
[correlationTS2, correlationTS_p2, correlationTS_rlo2, correlationTS_rup2]=movingCorrelation_c(Ma_corr2,sliding_w,1, 0.1);    


h7=plot(stacked_record_annual_Ma((80:110+add_yr),1), correlationTS2(:,1),'.-','LineWidth',4,'MarkerSize',14,'Color',[0.0,0.2,0.6]);  
h8=plot(stacked_record_annual_Ma((80:110+add_yr),1), correlationTS2(:,2),'--','LineWidth',4,'MarkerSize',14,'Color',[0.2,0.8,0.2]);  
               set(gca,...
          'linewidth',3,...
          'FontWeight','bold' );
    
   %%%%% Legend %%%%%%%%%%%%
    hs61='PC1 SAT';
    hs89='PC2 SAT';
    hs62='PC3 SAT';
    hs7='RICE SAT'; 
    hs8='EA SAT'; 
        
    hs7_c=['r (dD, ',hs7,')'];
    hs8_c=['r (dD, ',hs8,')'];

  h_leg=legend([ h7 h8 ],  hs7_c, hs8_c);
  set(h_leg, 'location', 'SouthWest','EdgeColor',[0 0 0],'FontSize',18,'FontWeight','bold' ); % ,'color','none'
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % letter
  letter_pos=[50 380 60 60];
  letter='c';
          TextBox = uicontrol('style','text');
          set(TextBox,'String',letter,'position',letter_pos,'FontWeight','bold','FontSize',34 ); % x position y position xsize ysize
          set(TextBox,'foregroundcolor', [0 0 0], ...
         'backgroundcolor', [1 1 1]);
 %%%%%%%%%%%%%%%%%%%%%%%%%
        set(gca,'XLim',[1977 2012])
        hl=hline([-0.5,0,0.5]);
        set(hl,'LineWidth',3);
        
%        xlabel('Year','FontWeight','bold','FontSize',18 );
         ylabel('Correlation','FontWeight','bold','FontSize',20 );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
% save fig.         
filename=['SAT_RICE_ncep_ncar_w_',num2str(sliding_w),'_', iso_str]; 
filedir ='C:\PHD\matlab_storage_of_output_files\figures\';
savefilename_c=strcat(filedir,filename);
% save as png 
orient landscape
% increase quality of picture by increase to -r300
export_fig('-png','-nocrop','-painters', '-depsc','-opengl', '-r170',savefilename_c); % PNG-nocrop' 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%% SAT PCs dD corr
% Fig. 6d
  
pol_nr31=1; %(polarity from Regression era)
pol_nr32=1;
pol_nr33=1;

fac_nr31=-1; % factor for comparison
fac_nr32=1;
fac_nr33=1; 

Ma_corr=NaN(le,6);
Ma_corr(:,1)=detrend(stacked_record_annual_Ma((80:110+add_yr),2));   
Ma_corr(:,2)=detrend(MA_PCs_save_SAT((1:31+add_yr),2))*pol_nr31*fac_nr31; % SAT PC1 SAM related    
Ma_corr(:,3)=detrend(MA_PCs_save_SAT((1:31+add_yr),3))*pol_nr32*fac_nr32;% SAT PC 2 PSA1 related  
Ma_corr(:,4)=detrend(MA_PCs_save_SAT((1:31+add_yr),4))*pol_nr33*fac_nr33;% SAT PC 3 PSA2 related
%edit movingCorrelation
[correlationTS, correlationTS_p, correlationTS_rlo, correlationTS_rup]=movingCorrelation_c(Ma_corr,sliding_w,1, 0.1);
    
 h=fig('units','inches','width',14,'height',5,'font','Helvetica','fontsize',16,'border','on');
   box on
   hold on   
 set(gca,'YLim',[-1 1])
            
 addpath C:\PHD\matlab_lib\Library
%   edit IPO_vline
    IPO_vline(15) 
   
       h31=plot(stacked_record_annual_Ma((80:110+add_yr),1), correlationTS(:,1),'.-','LineWidth',4,'MarkerSize',14,'Color',[1.0,0.0,1.0]);    % SAT PC1
       h32=plot(stacked_record_annual_Ma((80:110+add_yr),1), correlationTS(:,2),'.-','LineWidth',4,'MarkerSize',14,'Color',[0.1,0.1,0.1]);   % SAT PC2 
       h33=plot(stacked_record_annual_Ma((80:110+add_yr),1), correlationTS(:,3),'.-','LineWidth',4,'MarkerSize',14,'Color',[0.6,0.8,0.2]);     % SAT PC3
        set(gca,...
          'linewidth',3,...
          'FontWeight','bold' );
%%%%%%%%%%%%%%%
% legend
    hs31_c=['r (dD, ',hs61,' (',num2str(fac_nr31),'))'];
    hs32_c=['r (dD, ',hs89,')'];
    hs33_c=['r (dD, ',hs62,')'];
    h_leg=legend([ h31 h32  h33],  hs31_c, hs32_c , hs33_c); 
    set(h_leg, 'location', 'SouthWest','EdgeColor',[0 0 0],'FontSize',18,'FontWeight','bold' ); % ,'color','none'
      
%%%%%%%%%%%%%%%%%%%%
% legend
  letter_pos=[50 380 60 60];
  letter='d';
          TextBox = uicontrol('style','text');
          set(TextBox,'String',letter,'position',letter_pos,'FontWeight','bold','FontSize',34 ); % x position y position xsize ysize
          set(TextBox,'foregroundcolor', [0 0 0], ...
         'backgroundcolor', [1 1 1]);
 %%%%%%%%%%%%%%%%%%%%%%%%%
        set(gca,'XLim',[1977 2012])
        hl=hline([-0.5,0,0.5]);
        set(hl,'LineWidth',3);
        xlabel('Year','FontWeight','bold','FontSize',18 );
        ylabel('Correlation','FontWeight','bold','FontSize',20 );
%%%%%%%%%%%%
% save fig
      filename=['SAT_RICE_PCs_corr_w_',num2str(sliding_w),'_', iso_str];
      filedir ='C:\PHD\matlab_storage_of_output_files\figures\';
      savefilename_c=strcat(filedir,filename);
% save as png 
orient landscape
% increase quality of picture by increase to -r300
  export_fig('-png','-nocrop','-painters', '-depsc','-opengl', '-r170',savefilename_c); % PNG-nocrop'   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% SIC
  % (Fig 8b)
  % AS Ross Sea 
       close all
  % box based on correlation pattern 
  %  ~120W-180W, 65-75S
  
  HadISST_time=[1870:2014]'; % 145 
  load('C:\PHD\matlab_storage_of_output_files/Ma_dD_HadISST_SIC_box7_c.mat');
  load('C:\PHD\matlab_storage_of_output_files/Ma_dD_HadISST_SIC_box8.mat'); % eastern Ross Sea
  load('C:\PHD\matlab_storage_of_output_files/Ma_ADP_index_HadISST_SIC.mat'); % ADP index   ADP_index.m
        
   h=fig('units','inches','width',14,'height',5,'font','Helvetica','fontsize',16,'border','on');
   hold on      
   plot(MA_ADP_save(:,1),anomaly(MA_ADP_save(:,2))*(-1),'.-','LineWidth',4,'MarkerSize',14,'Color',[0.0,0.0,0.0]);
    
          box on
          set(gca,...
          'linewidth',3,...
          'FontWeight','bold' );
         
           set(gca,'YLim',[-3 3])
      %%%%%%%%%%%%%%            
      addpath C:\PHD\matlab_lib\Library
      IPO_vline(10)
       
      h1=plot(MA_ADP_save(:,1), anomaly(MA_ADP_save(:,2))*(-1),'*-k','LineWidth',4,'MarkerSize',14);
      h11=plot(HadISST_time(104:145), anomaly(HadISST_M_ice_box8(104:145))*(-1),'.-','LineWidth',4,'MarkerSize',14,'Color',[0.4,0.8,0.5]);      

       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  RICE iso
% annual
h15=plot(stacked_record_annual_Ma((1:110+add_yr),1),...
    anomaly(stacked_record_annual_Ma((1:110+add_yr),2)),'*-','Color',[0.5,0.5,0.5],'LineWidth',4,'MarkerSize',14);
           set(gca,...
          'linewidth',3,...
          'FontWeight','bold' );     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% legend
     hs1='ADP annual (-1)';
     hs11='E RS annual (-1)';
     hs15=[' {\delta}D annual'];
   
     h_leg=legend([ h1 h11 h15 ],  hs1, hs11 , hs15 );                     
     set(h_leg,'Position',[0.27 0.77 0.1 0.10] ,'EdgeColor',[0 0 0],'FontSize',18,'FontWeight','bold'); % ,'color','none'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% label
title_str='SIC';
  letter_pos=[1125 446 100 40];
          TextBox = uicontrol('style','text');
          set(TextBox,'String',title_str,'position',letter_pos,'FontWeight','bold','FontSize',26 ); % x position y position xsize ysize
          set(TextBox,'foregroundcolor', [0 0 0], ...
         'backgroundcolor', [1 1 1]); 
 %%%%%%%%%%%%%%%%%%%%%%%%%     
      set(gca,'XLim',[1977 2012])
      ylabel('Anomaly','FontWeight','bold','FontSize',20 );
            hl=hline(0.0);
            set(hl,'LineWidth',3);    
 %%%%%%%%%%%%%%%%%%%%%%%%%%
 % letter
       letter_pos=[50 380 60 60];  
       letter='b';
          TextBox = uicontrol('style','text');
          set(TextBox,'String',letter,'position',letter_pos,'FontWeight','bold','FontSize',34 ); % x position y position xsize ysize
          set(TextBox,'foregroundcolor', [0 0 0], ...
         'backgroundcolor', [1 1 1]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
%%% Save Fig. 
     filename=['SIC_RICE_HadISST_anomaly','_', iso_str];
     filedir ='C:\PHD\matlab_storage_of_output_files\figures\';
     savefilename_c=strcat(filedir,filename);
 % save as png 
     orient landscape
% increase quality of picture by increase to -r300
    export_fig('-png','-nocrop','-painters', '-depsc','-opengl', '-r170',savefilename_c); % PNG-nocrop'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%  %%%%%%%%%%%% SIC moving corr 
% Fig 8c
% SIC Regions 
           close all   
h=fig('units','inches','width',14,'height',5,'font','Helvetica','fontsize',16,'border','on');
hold on     

 load('C:\PHD\matlab_storage_of_output_files\HadISST_PCs_SIC_lim-150--30_-64_-75_1979-2014_annual.mat'); % use 
 MA_PCs_save_SIC=MA_PCs_save;
 MA_PCs_save_monthly_SIC=MA_PCs_save_monthly;
    
 %%%%
 add_yr=2; % 2- end 2011 ; 0- end 2009
 %%%%%%%%%%%%%%%%%%%%%%% RICE
le=length(stacked_record_annual_Ma(80:110+add_yr));
Ma_corr3=NaN(le,7);
% iso
Ma_corr3(:,1)=detrend(stacked_record_annual_Ma((80:110+ add_yr),2));
Ma_corr3((1:31+add_yr),2)=detrend(MA_ADP_save((8:38+add_yr),2)); % 1979=- % ADP index
Ma_corr3((1:31+add_yr),3)=detrend(MA_PCs_save_SIC((1:31+ add_yr),2)); %PC1, 1979
Ma_corr3((1:31+add_yr),4)=detrend(MA_PCs_save_SIC((1:31+ add_yr),3)); %PC2
Ma_corr3((1:31+add_yr),5)=detrend(HadISST_M_ice_box8(110:140+ add_yr)); % 1979 Eastern Ross Sea, starts 1870
Ma_corr3((1:31+add_yr),6)=detrend(MA_PCs_save_SIC((1:31+ add_yr),4)); % PC3
Ma_corr3((1:31+add_yr),7)=detrend(HadISST_M_ice_box7(110:140+ add_yr));
      
[correlationTS3, correlationTS_p3, correlationTS_rlo3, correlationTS_rup3]=movingCorrelation_c(Ma_corr3,sliding_w,1, 0.1);    

    plot([1979:2011], correlationTS3(:,1),'.-','LineWidth',4,'MarkerSize',14,'Color',[0.0,0.0,0.0]);     %   =-
    hold on
    box on
          set(gca,...
          'linewidth',3,...
          'FontWeight','bold' );
    
         set(gca,'YLim',[-1 1])
      %%%%%%%%%%%%%%
            shade_nr=1;
        if shade_nr==1
            addpath C:\PHD\matlab_lib\Library
            IPO_vline(10) % transitions as in England 2014
            % edit IPO_vline
        end
     
      h8=plot([1979:2011], correlationTS3(:,1),'.-','LineWidth',4,'MarkerSize',14,'Color',[0.0,0.0,0.0]);     %ADP
      h11=plot([1979:2011], correlationTS3(:,4),'-','LineWidth',4,'MarkerSize',14,'Color',[0.4,0.8,0.5]);     % 
      h14=plot([1979:2011], correlationTS3(:,6),'-','LineWidth',4,'MarkerSize',14,'Color',[0.1,0.1,0.7]);     %
       
      %%%%%%%%%%%%%%%%%%%%%% 
%            set(gca,'XLim',[1970 2010]) 
             set(gca,'XLim',[1977 2012])
            hl=hline([-0.5, 0, 0.5]);
            set(hl,'LineWidth',3);
  %       xlabel('Year','FontWeight','bold','FontSize',18 );
          ylabel('Correlation','FontWeight','bold','FontSize',20 );
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
      pol_nr21=1; % Polarity defined in Fig 2 and 7 
      pol_nr22=-1;
      pol_nr23=-1;
      
      fac_nr21=-1; % switch for comparison
      fac_nr22=1;
      fac_nr23=1;              
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % legend
            hs8='ADP SIC';
            hs11='E RS SIC';
            hs14='AS/RS SIC';
                   
            hs9='PC1 SIC';
            hs10='PC2 SIC';
            hs12='PC3 SIC';
            
               hs8_c=['r (dD, ',hs8,')'];
               hs9_c=['r (dD, ',hs9,'(',num2str(fac_nr21),'))'];
               hs10_c=['r (dD, ',hs10,')'];
               hs11_c=['r (dD, ',hs11,')'];
               hs12_c=['r (dD, ',hs12,')'];
               hs14_c=['r (dD, ',hs14,')'];   
            
    h_leg=legend([ h8 h11 h14 ],  hs8_c, hs11_c, hs14_c );    
    set(h_leg, 'location', 'NorthEast','EdgeColor',[0 0 0],'FontSize',18,'FontWeight','bold'); % ,'color','none'        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% letter           
             letter='c';
             letter_pos=[50 380 60 60];  
 
          TextBox = uicontrol('style','text');
          set(TextBox,'String',letter,'position',letter_pos,'FontWeight','bold','FontSize',34 ); % x position y position xsize ysize
          set(TextBox,'foregroundcolor', [0 0 0], ...
         'backgroundcolor', [1 1 1]); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save fig.
      filename=['SIC_RICE_HadISST_corr_moving_w_',num2str(sliding_w),'_', iso_str];      
      filedir ='C:\PHD\matlab_storage_of_output_files\figures\';
      savefilename_c=strcat(filedir,filename);

 % save as png 
orient landscape
% increase quality of picture by increase to -r300
export_fig('-png','-nocrop','-painters', '-depsc','-opengl', '-r170',savefilename_c); % PNG-nocrop' 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% moving corr SIC
% Fig 8d
% PC1, PC2 and PC3 
           close all
        
   h=fig('units','inches','width',14,'height',5,'font','Helvetica','fontsize',16,'border','on');
   hold on 
   box on
          set(gca,...
          'linewidth',3,...
          'FontWeight','bold' );
    
         set(gca,'YLim',[-1 1])
      %%%%%%%%%%%%%%
        shade_nr=1;
        if shade_nr==1
            addpath C:\PHD\matlab_lib\Library
            IPO_vline(15)
            hold on
        end
      
      h9=plot([1979:2011], correlationTS3(:,2)*pol_nr21*fac_nr21,'.-','LineWidth',4,'MarkerSize',14,'Color',[0.6,0.2,0.2]);     % PC1
      h10=plot([1979:2011], correlationTS3(:,3)*pol_nr22*fac_nr22,'.-','LineWidth',4,'MarkerSize',14,'Color',[0.0,0.6,0.1]);     % PC2
      h12=plot([1979:2011], correlationTS3(:,5)*pol_nr23*fac_nr23,'--','LineWidth',4,'MarkerSize',14,'Color',[1.0,0.1,1.0]);     % PC3

          set(gca,'XLim',[1977 2012])
          hl=hline([-0.5, 0, 0.5]);
          set(hl,'LineWidth',3);
          xlabel('Year','FontWeight','bold','FontSize',18 );
          ylabel('Correlation','FontWeight','bold','FontSize',20 );       
%%%%%%%%%
% legend
        h_leg=legend([  h9 h10 h12 ],  hs9_c, hs10_c , hs12_c );    
        set(h_leg, 'location', 'SouthWest','EdgeColor',[0 0 0],'FontSize',18,'FontWeight','bold'); % ,'color','none'
%%%%%%%%%    
% letter            
             letter='d';
             letter_pos=[50 380 60 60];  
 
          TextBox = uicontrol('style','text');
          set(TextBox,'String',letter,'position',letter_pos,'FontWeight','bold','FontSize',34 ); % x position y position xsize ysize
          set(TextBox,'foregroundcolor', [0 0 0], ...
         'backgroundcolor', [1 1 1]); 
%%%%%%%%%%%%%%%%%%%          
% save fig.        
      filename=['SIC_PCs_RICE_HadISST_corr_moving_w_',num2str(sliding_w),'_', iso_str]; 
      filedir ='C:\PHD\matlab_storage_of_output_files\figures\';
      savefilename_c=strcat(filedir,filename);

% save as png 
orient landscape
% increase quality of picture by increase to -r300
export_fig('-png','-nocrop','-painters', '-depsc','-opengl', '-r170',savefilename_c); % PNG-nocrop'   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% z500
% Fig 5b 
load('C:\PHD\matlab_storage_of_output_files\ERA-Interim_PCs_z500_lim0-360_-20_-90_1979-2014_annual.mat') %z500
MA_PCs_save_z500=MA_PCs_save;  
MA_PCs_save_monthly_z500=MA_PCs_save_monthly;
clear MA_PCs_save MA_PCs_save_monthly  
%%%%%%%%%%%%%%%%%%%%%%% RICE
pol_nr11=1; % polarity positive SAM, as displayed in Fig 2
pol_nr12=1; % PSA patterns as in Kidson 1988 his fig 4b,c 
pol_nr13=-1;

fac_nr11=-1; 
fac_nr12=1;
fac_nr13=1;

le=length(stacked_record_annual_Ma(80:110+add_yr)); % 1979-2011
Ma_corr=NaN(le,6);
Ma_corr4(:,1)=detrend(stacked_record_annual_Ma((80:110+add_yr),2));
Ma_corr4((1:31+add_yr),2)=detrend(MA_PCs_save_z500((1:31+add_yr),2))*pol_nr11*fac_nr11; % z500 PC1 SAM related  (factor from Regression era)
Ma_corr4((1:31+add_yr),3)=detrend(MA_PCs_save_z500((1:31+add_yr),3))*pol_nr12*fac_nr12; % z500 PC 3 PSA2 related  
Ma_corr4((1:31+add_yr),4)=detrend(MA_PCs_save_z500((1:31+add_yr),4))*pol_nr13*fac_nr13; % z500 PC 3 PSA1 related

%edit movingCorrelation

 [correlationTS4, correlationTS4_p, correlationTS4_rlo, correlationTS4_rup]=movingCorrelation_c(Ma_corr4,sliding_w,1, 0.1);
 
 h=fig('units','inches','width',14,'height',5,'font','Helvetica','fontsize',16,'border','on');
   box on
   hold on 
   set(gca,'YLim',[-1 1])
 %%%%%%%%%%%%%%

        if shade_nr==1
            
            addpath C:\PHD\matlab_lib\Library
            IPO_vline(10)
            hold on
        
        end
    h1=plot(stacked_record_annual_Ma((80:110+add_yr),1), correlationTS4(:,1),'.-','LineWidth',4,'MarkerSize',14,'Color',[0.0,0.2,0.6]);     %   
    h2=plot(stacked_record_annual_Ma((80:110+add_yr),1), correlationTS4(:,2),'-','LineWidth',4,'MarkerSize',14,'Color',[0.2,0.4,0.0]);     
    h3=plot(stacked_record_annual_Ma((80:110+add_yr),1), correlationTS4(:,3),'-','LineWidth',4,'MarkerSize',14,'Color',[0.8,0.1,0.0]);
    
               set(gca,...
          'linewidth',3,...
          'FontWeight','bold' ); 
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% legend
%     hs1='PC1 z500';
%     hs2='PC2 z500';
%     hs3='PC3 z500';
    
    hs1='SAM';
    hs2='PSA1';
    hs3='PSA2';
    
    hs1_c=['r (dD, ',hs1,' (',num2str(fac_nr11),'))'];
    hs2_c=['r (dD, ',hs2,')'];
    hs3_c=['r (dD, ',hs3,')'];
 
    
  h_leg=legend([ h1 h2 h3],  hs1_c, hs2_c, hs3_c); 
  set(h_leg, 'location', 'SouthWest','EdgeColor',[0 0 0],'FontSize',18,'FontWeight','bold' ); % ,'color','none'
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % letter
  letter_pos=[90 380 60 60];
  letter='b';
          TextBox = uicontrol('style','text');
          set(TextBox,'String',letter,'position',letter_pos,'FontWeight','bold','FontSize',34 ); % x position y position xsize ysize
          set(TextBox,'foregroundcolor', [0 0 0], ...
         'backgroundcolor', [1 1 1]);
 %%%%%%%%%%%%%%%%%%%%%%%%%
        set(gca,'XLim',[1977 2012])
        hl=hline([-0.5,0,0.5]);
        set(hl,'LineWidth',3);
        xlabel('Year','FontWeight','bold','FontSize',18 );
        ylabel('Correlation','FontWeight','bold','FontSize',20 );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save fig
filename=['z500_RICE_ERA_interim_mov_corr_w_',num2str(sliding_w),'_', iso_str];       
filedir ='C:\PHD\matlab_storage_of_output_files\figures\';
savefilename_c=strcat(filedir,filename);
% save as png 
orient landscape
% increase quality of picture by increase to -r300
export_fig('-png','-nocrop','-painters', '-depsc','-opengl', '-r170',savefilename_c); % PNG-nocrop' 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 %% Fig. 4 All PCs (z500, SAT, SIC)
 alt_nr=3; % (1) PC1 z500 SAT SIC (2) PC2 z500 and SIC, PC3 SAT (3)
   % (b) anomaly series 
 h=fig('units','inches','width',14,'height',5,'font','Helvetica','fontsize',16,'border','on');
   box on
   hold on   
%  plot(stacked_record_annual_Ma((49:110),1), correlationTS(:,1),'.-','LineWidth',4,'MarkerSize',14,'Color',[0.0,0.0,0.0]);  % not visable
     set(gca,'XLim',[1977 2012])
     set(gca,'YLim',[-1.5 1.5])
 %%%%%%%%%%%%%%
        if shade_nr==1
            addpath C:\PHD\matlab_lib\Library
            IPO_vline(10)
            hold on
        end

           set(gca,...
          'linewidth',3,...
          'FontWeight','bold' );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
 % fac_nr not used here, that is no sign is switched for comparision, thus all PCs are displayed with their positive polarity   
        
 % z500, all positive polarity
   if alt_nr==1
    pol_nr=1;
    fac_nr=1;
    h1=plot(MA_PCs_save_z500((1:33),1), MA_PCs_save_z500((1:33),2)*pol_nr*fac_nr,'.-','LineWidth',4,'MarkerSize',14,'Color',[0.0,0.2,0.6]); %z500 EOF1
    hs1_c='z500 PC1';
    
   elseif alt_nr==2
      
    pol_nr=1;
    fac_nr=1; 
    h1=plot(MA_PCs_save_z500((1:33),1), MA_PCs_save_z500((1:33),3)*pol_nr*fac_nr,'-','LineWidth',4,'MarkerSize',14,'Color',[0.0,0.2,0.6]);  %z500 EOF2
    hs1_c='z500 PC2';
    
   elseif alt_nr==3
     
    pol_nr=-1;
    fac_nr=1; 
    h1=plot(MA_PCs_save_z500((1:33),1), MA_PCs_save_z500((1:33),4)*pol_nr*fac_nr,'-','LineWidth',4,'MarkerSize',14,'Color',[0.0,0.2,0.6]);  %z500 EOF3
    hs1_c='z500 PC3';
    
   end
      
 % SAT   
   if alt_nr==1
    pol_nr=1;
    fac_nr=1;
    h2=plot(MA_PCs_save_z500((1:33),1), MA_PCs_save_SAT((1:33),2)*pol_nr*fac_nr,'-','LineWidth',4,'MarkerSize',14,'Color',[0.8,0.1,0.0]); % SAT EOF1
    hs2_c='SAT PC1';
   elseif alt_nr==2
    pol_nr=1;
    fac_nr=1;
    h2=plot(MA_PCs_save_z500((1:33),1), MA_PCs_save_SAT((1:33),4)*pol_nr*fac_nr,'-','LineWidth',4,'MarkerSize',14,'Color',[0.8,0.1,0.0]); % SAT EOF3
    hs2_c='SAT PC3';
   elseif alt_nr==3
    pol_nr=1;
    fac_nr=1;
    h2=plot(MA_PCs_save_z500((1:33),1), MA_PCs_save_SAT((1:33),3)*pol_nr*fac_nr,'-','LineWidth',4,'MarkerSize',14,'Color',[0.8,0.1,0.0]); % SAT EOF2
    hs2_c='SAT PC2';
   end

   % SIC
   if alt_nr==1
    pol_nr=1;
    fac_nr=1;
    h3=plot(MA_PCs_save_z500((1:33),1), MA_PCs_save_SIC((1:33),2)*pol_nr*fac_nr,'--','LineWidth',4,'MarkerSize',14,'Color',[0.1,0.1,0.1]);  %SIC EOF1
    hs3_c='SIC PC1';
   elseif alt_nr==2
    pol_nr=-1;
    fac_nr=1;
    h3=plot(MA_PCs_save_z500((1:33),1), MA_PCs_save_SIC((1:33),3)*pol_nr*fac_nr,'--','LineWidth',4,'MarkerSize',14,'Color',[0.1,0.1,0.1]);  %SIC EOF2
    hs3_c='SIC PC2';
   elseif alt_nr==3
    pol_nr=-1;
    fac_nr=1;
    h3=plot(MA_PCs_save_z500((1:33),1), MA_PCs_save_SIC((1:33),4)*pol_nr*fac_nr,'--','LineWidth',4,'MarkerSize',14,'Color',[0.1,0.1,0.1]);  %SIC EOF3
    hs3_c='SIC PC3';
   end

  h_leg=legend([ h1 h2 h3],  hs1_c, hs2_c, hs3_c); 
  set(h_leg, 'location', 'NorthWest','EdgeColor',[0 0 0],'FontSize',18,'FontWeight','bold' ); % ,'color','none'      
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % letter
  letter_pos=[70 380 60 60];
  if alt_nr==1
    letter='a';
  elseif alt_nr==2
    letter='b'; % xlabel('Year','FontWeight','bold','FontSize',18 );
  elseif alt_nr==3
    letter='c'; xlabel('Year','FontWeight','bold','FontSize',18 );  
  end
          TextBox = uicontrol('style','text');
          set(TextBox,'String',letter,'position',letter_pos,'FontWeight','bold','FontSize',34 ); % x position y position xsize ysize
          set(TextBox,'foregroundcolor', [0 0 0], ...
         'backgroundcolor', [1 1 1]);
%%%%%%%%%%%%%%%%%%%%%%%%%
% save fig  
filename=['z500_SAT_SIC_PC_ts', iso_str,'_alt_',num2str(alt_nr)]; 
filedir ='C:\PHD\matlab_storage_of_output_files\figures\';
savefilename_c=strcat(filedir,filename);
 % save as png 
orient landscape
% increase quality of picture by increase to -r300
export_fig('-png','-nocrop','-painters', '-depsc','-opengl', '-r170',savefilename_c); % PNG-nocrop' 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Regression_spatial_HadISST_SIC_DE18_c4.m
% Matlab 2017a
% Daniel Emanuelsson 2018
% SIC Regression
% Github version 1   [uploaded]
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subfunctions
% *fig.m fileexchange
% *export_fig fileexchange
%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
% type_nr=2;    % (1) annual (2) monthly

param_nr=2;% (1-SST/ 2 SIC)

name='SIC';
site='RICE';
 
 set_yr=1;
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if set_yr==1    
        yr_s=1979;
        yr_e=2014;     
   elseif set_yr==8       
        yr_s=2000;
        yr_e=2013;
    end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
 show_colorbar=1;
 sea_nr=1;
 season='annual';

 
 %%%%%%%%%%%%%%%%%%%%
 iso_nr=15;  
 %%%%%%%%%%%%%%%%%%% 
if iso_nr==7 %***************
    iso='PCs';
elseif iso_nr==9 %***************
    iso=[]; % just trend
elseif iso_nr==14 %***************
    iso='CP EP index';
elseif iso_nr==15 %***************
    iso='SIC PCs';
end
 
show_maximum_point=0;
 
show_title=1; % (1/0)
area_2_box=0 ; % 0/1 for RICE SST corr Box and text for Area 2
  
% if  strcmp(name,'SIC')==1
  
  lat1=-90;
  lat2=-50;
  box_use=0;  % (0/1) 

  lon1=-300;
  lon2= 60;
  
  lock_scalebar=1;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SST_dataset=1; % (1) HadISST
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HadISST data
% version 1.1 Rayner et al. 2003


addpath C:\PHD\HadISST_data\ % Ienovo
%    ncdisp('HadISST_ice_c.nc');
%    ncdisp('HadISST_sst.nc');

 HadISST_time=ncread('HadISST_ice_c.nc','time'); %units         = 'days since 1870-1-1 0:0:0'
 % missing values  missing_value = -1e+30 ( looks more like -1000)\
 %monthly values
 
 HadISST_time_c=HadISST_time+datenum(1870,1,1);  % to 2015 month 7
 
 date_vec=datevec(double(HadISST_time_c)); 
 yyyy=date_vec(:,1); 
 mm=date_vec(:,2);
 HadISST_year_num=yyyy+(mm/12-1/12); 

 HadISST_time_bnds =ncread('HadISST_ice_c.nc','time_bnds');
 HadISST_lat=ncread('HadISST_ice_c.nc','latitude');
 HadISST_lon=ncread('HadISST_ice_c.nc','longitude');
 

 if strcmp(name,'SIC')==1  
%   name='ice'; 
   HadISST_ice=ncread('HadISST_ice_c.nc','sic'); % 360x180x1728 lon x lat x time
   M=HadISST_ice;
 
 end
   
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % 4. 0 Seasonal and Annual
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %
 %
 % HadISST_count=1705; % that overlaps with RICE record 1979-2011
HadISST_count=find(HadISST_year_num==yr_e)+11; % values in months
% HadISST_start=289; %1894
HadISST_start=find(HadISST_year_num==yr_s); % 1950 hardly any data before this
% HadISST_start=1309; %1979 satelite
% (HadISST_count-HadISST_start)/12

% era_count=1693; 
% div=118;
div=ceil((HadISST_count-HadISST_start)/12); %1979

%%%%%%%%%%%%%%%%%%%%%%%

% mm_in=mm(HadISST_start:HadISST_count);  % seasonal index

mm_in_c=[HadISST_start:HadISST_count]';
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % cosweight
% S_lon=17;S_lat=170;%check
% M ( S_lon, S_lat-5,2) 
%  M=cosweight_c(M,HadISST_lat); % UoW function
%  M ( S_lon, S_lat-5,2)


HadISST_year_num(HadISST_start);
HadISST_year_num(HadISST_count);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% 
     
if iso_nr==7

load('C:\PHD\matlab_storage_of_output_files\ERA-Interim_PCs_z500_lim0-360_-20_-90_1979-2014_annual'); % 2014

            time_c=MA_PCs_save_monthly(:,1);
            yr_2=find( time_c==yr_e)+11;
            yr_1=find( time_c==yr_s); 
               
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            PC_nr=2; % (2) SAM (3) PSA1 (4) PSA2
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
             label_2='';
             iso='z500';
   
            if PC_nr==2 %SAM
                    letter='d';
                    fact=1; % to corespond to SIC PC regression fig. 2014   
            elseif PC_nr==3  %PSA1+
                    letter='e';
                    fact=1;   %PSA1   % to corespond to SIC PC regression fig. 2014     
            elseif PC_nr==4
                    letter='f';
                    fact=-1; % to corespond to SIC PC regression fig.2014
            end
                      
 y=MA_PCs_save_monthly((yr_1:yr_2),PC_nr)*fact;
              
  elseif iso_nr==9
      
      
            wr_yr=[1979:2014]';
            yr_1=find(wr_yr==yr_s);       
            yr_2=find(wr_yr==yr_e);
            
            yc=[1:HadISST_count]';
            y=yc(mm_in_c);
                             
            PC_nr=999; % just to get a number      
            label_2='yr^-^1';
            

                     if set_yr==1
                         letter='a';
                     elseif set_yr==8
                         letter='a'; 
                     end
                 
                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
    elseif iso_nr==14      %CP EP index             
       
        el_nino_type=1; %  (1) CP (2) EP 
                
       if el_nino_type==1
           
%           load('C:\PHD\matlab_storage_of_output_files\HadISST_PCs_SST_CP_lim-70_120_20_-20_1950-2014_annual_c.mat') % 1950
          load('C:\PHD\matlab_storage_of_output_files\HadISST_PCs_SST_CP_lim-70_120_20_-20_1900-2014_annual_c.mat') % 1900  fact (-1) used
                   
       elseif el_nino_type==2
           
           %load('C:\PHD\matlab_storage_of_output_files\HadISST_PCs_SST_EP_lim-70_120_20_-20_1950-2014_annual_c.mat')
           load('C:\PHD\matlab_storage_of_output_files\HadISST_PCs_SST_EP_lim-70_120_20_-20_1900-2014_annual_c.mat')      
       end
       
                   time_c=MA_PCs_save_monthly(:,1);
                    yr_2=find( time_c==yr_e)+11;
                    yr_1=find( time_c==yr_s); 
       
                PC_nr=2;% leading PC
                fact=-1; % 1950 +1, 1900 -1
                y=MA_PCs_save_monthly((yr_1:yr_2),PC_nr)*fact;
                PC_nr=999; % just to get a number for enso    
                letter='b';      
                label_2='°C s.d.^-^1';
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    elseif iso_nr==15   % SIC PCs
        
        SIC_PCs_alt=5;%%%%%%%%%%%%%%%
        
          load('C:\PHD\matlab_storage_of_output_files\HadISST_PCs_SIC_lim-150--30_-64_-75_1979-2014_annual.mat'); % use 
          
  
                    time_c=MA_PCs_save_monthly(:,1);
                    yr_2=find( time_c==yr_e)+11;
                    yr_1=find( time_c==yr_s); 
                %%%%%%%%%%%%%%%%%%%%%%    
                PC_nr=2;%  PC# (2) PC1 Leading SAM-related (3) PC2 PSA1-related (4) PC3
                %%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%
                 if PC_nr==2
                    letter='a'; 
                    fact=1; %%%%%%%%%%%%% leave SIC PCs unchanged                 
                 elseif PC_nr==3
                    letter='b';
                    fact=-1; %%%%%%%%%%%%%
                 elseif PC_nr==4
                    letter='c'; 
                    fact=-1; %                 
                 end
                 
                y=MA_PCs_save_monthly((yr_1:yr_2),PC_nr)*fact;
                 
                label_2='SIC s.d.^-^1';
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                  
           
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% regression 
% regstats
 
%         if (set_yr==1 || set_yr==8  || set_yr==9) 
        load('C:\PHD\matlab_storage_of_output_files\HadISST_clim_sic_1979_2014.mat'); % clim removed saaved from eof_code_run_annual_seasonal_c17.m file
        ntim=432;
%         end

    
        nlat=180;
        nlon=360;
        
    
    HadISST_sst_c5 = reshape(HadISST_sst_c4,ntim, nlat, nlon); 
    HadISST_sst_c6=permute(HadISST_sst_c5,[3 2 1]);
  %  X=HadISST_sst_c6(:,:,(yr_1:yr_2));
  
            
            if set_yr==1 
                X=HadISST_sst_c6(:,:,(1:end)); 
            elseif set_yr==8
                X=HadISST_sst_c6(:,:,(253:end-12)); 
            end
    
    
[A B C]=size(X);
% b=zeros(2, A,B);bINT=zeros(2,2,A,B);r=zeros(324,A,B);rINT=zeros(324,2,A,B);
b=zeros(2, A,B);
bINT=zeros(2,2,A,B);
r=zeros(C,A,B); % residual
rINT=zeros(C,2,A,B);
STATS=zeros(1,4,A,B);  % 3 p-value

 p_all=zeros(A,B,1);  s_all=zeros(A,B,1);  % just one layer needed?
 
%  s_all_c_in=s_all;


for i=1:A
    for j=1:B

       
            
%%%%%%%%%%%%%%%%%
% Replace NaNs
            if  isnan(X(i,j,:))==1
               X(i,j,:)=1000;
            end
%%%%%%%%%%%%%%%%%%%

            
            s=regstats(squeeze(X(i,j,:)),y,'linear','all');
            s.tstat;
            s_all(i,j,:)=s.tstat.beta(2); % beta- regresion coeficient 
            p_iso_rec=s.tstat.pval;
            p_all(i,j,:)= p_iso_rec(2); % The second element is the p-value for a significantly non-zero slope.
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%% Put (z500, u850) field first
            %%%%% check slop in at the end of file

    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Figure
%%%%%%%%%%%%%%%%%%%%%%%
%%% Monthly values

proj_nr=1;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
proj= 'stereo';
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if   iso_nr==9 % TREND  per decade  % monthly                       
%                SIC in % /decade
        s_all_c1=s_all*10*12*100;
    else    % unit per standard dev.                   
        s_all_c1=s_all;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
[c_max, c_min, c_limit]=Ma_max_min(s_all_c1);

c_max_c=c_max;

c_limit_c=c_limit;

    if lock_scalebar==1      
        if  iso_nr==7 || iso_nr==15 % PCs
           c_limit=0.1;
           adj_nr=25;
        elseif  iso_nr==9   % trend SIC
           c_limit=0.15*100;
           adj_nr=25;        
        end        
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% Significance level    
s_all_c=s_all_c1;
if  iso_nr==9    % trend
    p_level_rc=1; 
else
    p_level_rc=0.05;   
end
%p_level_rc=0.01;

s_all_c(p_all>=p_level_rc)=NaN; 

%%%%%%%%%%%%%
[wa, wb]=size(s_all_c);

fig('units','inches','width',14,'height',9,'font','Helvetica','fontsize',16,'border','on');

if   iso_nr==9  
    
        lat2=-50; 
        adj_nr=9.5;
  
elseif  iso_nr==14  || iso_nr==15

        lat1=-50;    
        lat2=-90; 
        adj_nr=0;
  
    
end

% Define map axes and set map properties
      axesm( proj,'MapLatLimit',[lat1 lat2],'Grid','on','ParallelLabel','on','Frame','on',... %
       'MeridianLabel','on','FontWeight','bold','FontSize',18,...
       'mlabellocation',[0:30:180,0:-30:-180]); 
      
  
    set(gca,'box','off','XColor',[1,1,1],'YColor',[1,1,1]);  % turns off border box and axes
    wr_cc=squeeze(s_all_c(:,:,1))';
    hSurf=surfm(double(HadISST_lat),double(HadISST_lon),wr_cc);

    hold on
    colormap(b2r(-c_limit,c_limit));

%%%%%%%%
% change colorbar
br=colorbar;

%the position arguments are [xposition yposition width height].

        pos=get(br,'Position');
        pos(1)=pos(1)-0.0;      %%%%%%%%%%%%%%%%%%%% x-position of colorbar   (+) --->>   +0.02;
                       
            if iso_nr==9  ||  iso_nr==14               
                pos(1)=pos(1)+0.005; 
                pos(2)=pos(2)+0.07; 
                pos(3)=pos(3)- 0.0;  % widthcolorbar
                pos(4)=pos(4)-0.15;  % height colorbar

            elseif  iso_nr==7 ||  iso_nr==15 
                pos(1)=pos(1)+0.014; 
                pos(2)=pos(2)+0.170;
                pos(3)=pos(3)- 0.0;  % widthcolorbar
                pos(4)=pos(4)-0.35;  % height colorbar
                
            elseif param_nr==2
                pos(2)=pos(2)+0.045;   
                pos(3)=pos(3)- 0.0;  % widthcolorbar
                pos(4)=pos(4)-0.15;  % height colorbar   
            end

        
        set(br,'Position',pos)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Color of shading
% color_alt=2;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if iso_nr==7 || iso_nr==15 || iso_nr==14
    color_alt=4; 
elseif  iso_nr==9  
    color_alt=1; 
end

if color_alt==1
    colormap(brewermap(256,'*RdBu'));
elseif color_alt==2
    colormap(flipud(cbrewer('div','Spectral',10))); 
elseif color_alt==3
    colormap(flipud(cbrewer('div','PiYG',20)));    
elseif color_alt==4    
    colormap(cbrewer('div','BrBG',20));    
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Significance level contour

    if iso_nr==7 || iso_nr==9  || iso_nr==14 ||  iso_nr==15
        co_pa=[.1 .1 .1];
        
        line_w_nr=2;
        
    else
        co_pa=[.9 .9 .9];
        line_w_nr=2;
    end

    h1= contourm( double(HadISST_lat),double(HadISST_lon),squeeze(p_all(:,:,1))', [0.05],'--','ShowText','off',...
    'Linecolor',co_pa,'LineWidth', line_w_nr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coastline
coast_nr=0; %(1/0) on/off

if coast_nr==1
    load coast
    %%%%%%%%%%%%%
    % to be able to use coastline for continents in combination with ant
    % coastline and grounding line from bedmap
    in_c=find(lat<-60);
    lat_cr=lat;
    lon_cr=long;
    lat_cr(in_c)=NaN;
    lon_cr(in_c)=NaN;


%      color_code=[.1 .1 .1];
       color_code=[.0 .0 .0];  % black      
%      color_code=[.6 .6 .6];  % gray


    %%%%%%%%%%%%%
    plot3m(lat_cr,lon_cr,'-','LineWidth', 2, 'color',color_code);
    %plot3m(lat(in_c),long(in_c),'.k','MarkerSize', 1);
    %%%%%%%%%%%%%

    addpath C:\PHD\matlab_mapping
    bedmap2('gl','LineWidth', 1.5, 'color',color_code);
    bedmap2('patchshelves','LineWidth', 1.0,'facecolor','none','frame','on','edgecolor',color_code);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
end

%%%%%%%%%%%%%%%%%%%%
% Lat lon labels
label_vert=-300;
if param_nr==2 && iso_nr==9    
    label_across=-55;
    label_vert=-298;
elseif param_nr==2
     %label_across=-85;  
     label_across=-55; 
end
    mlabel('MLabelParallel',label_across,'PLabelLocation',[-75 -60 -45 -30 -15  0 15 ],'PLabelMeridian',label_vert) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if  iso_nr==7
                 
         label_1='°C s.d.^-^1'; % per nomralized unit standard deviation
         
         
         if param_nr==1   
            x_lab=0.920;y_lab=0.7;
         elseif param_nr==2
            label_1='SIC';   
            %x_lab=1.04;y_lab=0.41;
            x_lab=1.09;y_lab=0.45;
         end

        if PC_nr==2 
         pc_str='SAM';

        elseif PC_nr==3 && iso_nr==7
            pc_str='PSA1';

        elseif PC_nr==4 && iso_nr==7
            
            if fact==1
            pc_str='PSA2';
            elseif fact==-1
            pc_str='PSA2';
            end
        end
                                
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
      elseif   iso_nr==9 
%             if  iso_nr==9
                pc_str='Trend';
                
                if param_nr==1
              	label_1='°C decade^-^1';
                x_lab=0.920;y_lab=0.74;
                
                elseif param_nr==2
                label_1='% decade^-^1';
                
                    alt_cp1=2;%%%%
                    if alt_cp1==1
                        x_lab=0.920;y_lab=0.64;
                    elseif alt_cp1==2
                        x_lab=1.075;y_lab=0.34;
                    elseif alt_cp1==3
                        x_lab=1.025;y_lab=0.10;
                    end
                    
                end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                    
      elseif  iso_nr==14
                            if el_nino_type==1
                                pc_str='CP';
                            elseif el_nino_type==2
                                pc_str='EP';
%                             elseif el_nino_type==12
%                                 pc_str='CP alt 2';
                            end
                            
                            
                               if param_nr==1
                                label_1='°C s.d.^-^1'; % per nomralized unit standard deviation     
                                x_lab=0.920;y_lab=0.67;  
                        
                              elseif param_nr==2
                                label_1='SIC';   
                                %x_lab=1.04;y_lab=0.42;
                                
                                x_lab=1.1;y_lab=0.44;
                                
                               end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                          
    elseif  iso_nr==15  
                        label_size=26;
                        if PC_nr==2
                        %pc_str='PC1';
                        pc_str='EOF1';
                                    if SIC_PCs_alt==1
                                        variance_exp=17;
                                    elseif SIC_PCs_alt==2
                                        variance_exp=17;
                                    elseif SIC_PCs_alt==5
                                        variance_exp=17;  
                                    elseif SIC_PCs_alt==4
                                        variance_exp=15;
                                    end
                                    
                                    a2=axestext_c(0.95, +0.96, [num2str(variance_exp),'%'] );
                                    set(a2,'FontWeight','bold','FontSize',label_size);
                                                                 
                                    % SIC Sector for EOF         
                                            lon_b_c=[-30 -40 -50 -60 -70 -80 -90 -100 -110 -120 -130 -140 -150   -150 -140 -130 -120 -110 -100 -90 -80 -70 -60 -50 -40 -30 -30 ];
                                            lat_b_c=[-75 -75 -75 -75 -75 -75 -75  -75  -75  -75  -75  -75  -75   -64   -64  -64  -64  -64  -64 -64 -64 -64 -64 -64 -64 -64 -75 ];
                                            plotm(lat_b_c,lon_b_c,'-','LineWidth',3,'Color',[1.0 0.0 .4])                      
                        
                        elseif PC_nr==3
                        %pc_str='PC2';
                        pc_str='EOF2';
                                     if SIC_PCs_alt==1
                                        variance_exp=11;
                                     elseif SIC_PCs_alt==2
                                        variance_exp=11; 
                                    elseif SIC_PCs_alt==5  %%%
                                        variance_exp=12;
                                    elseif SIC_PCs_alt==4
                                        variance_exp=10;
                                     end
                                    a2=axestext_c(0.95, +0.96, [num2str(variance_exp),'%'] );
                                    set(a2,'FontWeight','bold','FontSize',label_size);
                                    
                       elseif PC_nr==4
                        %pc_str='PC3';
                        pc_str='EOF3';
                                     if SIC_PCs_alt==1
                                        variance_exp=11;
                                     elseif SIC_PCs_alt==2
                                        variance_exp=8; 
                                     elseif SIC_PCs_alt==5 %%%%
                                        variance_exp=8;
                                     elseif SIC_PCs_alt==4
                                        variance_exp=10;
                                     end
                                    a2=axestext_c(0.95, +0.96, [num2str(variance_exp),'%'] );
                                    set(a2,'FontWeight','bold','FontSize',label_size);              
                        
                        
                        end
                               
                               if param_nr==2
%                                 label_1='SIC s.d.^-^1';   
                                label_1='SIC';   
%                                 x_lab=1.04; y_lab=0.41;
                               x_lab=1.09;y_lab=0.45;
                               end      
            
    end
    
%%%%%%%%%%%%%%%
% title
    season=[];
    type_str='monthly';

%      t4=title([ pc_str  ,' SIC regression ',season,' ',type_str,' ',num2str(yr_s),'-',num2str(yr_e) ]); 
     t4=title([ pc_str  ,' SIC regression ',type_str,' ',num2str(yr_s),'-',num2str(yr_e) ]);      
           

%%%%%%%%% Move title
        pos=get(t4,'Position');
        
        if  iso_nr==9  
%             pos(2)=pos(2)-0.24;
              pos(2)=pos(2)-0.06;
        elseif (iso_nr==7 || iso_nr==14  || iso_nr==15 ) && param_nr==2             
            pos(2)=pos(2)-0.06;

        else
            pos(2)=pos(2)-0.16;  
        end
        
        set(t4,'Position',pos,'FontSize',22, 'FontWeight', 'bold');    
%%%%%%%%%%%%%%%%%%%%
% letter
        if  iso_nr==9 && param_nr==2
            letter_pos= [350 685 50 50];
        elseif  ( iso_nr==7 || iso_nr==14|| iso_nr==15 ) && param_nr==2 
            sic_alt=2;%%%%%%%%%%%%%%%%%%%%%%
            if  sic_alt==1
                letter_pos= [270 545 50 50];
            elseif sic_alt==2
                letter_pos= [380 675 50 50];
            end
        else
            letter_pos= [270 560 50 50];   
        end

         letter_size=36;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%

          TextBox = uicontrol('style','text');
          set(TextBox,'String',letter,'position',letter_pos,'FontWeight','bold','FontSize',letter_size ); % x position y position xsize ysize
          set(TextBox,'foregroundcolor', [0 0 0], ...
         'backgroundcolor', [1 1 1]);
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % label  SIC
     if param_nr==2 && ( iso_nr==7 ||  iso_nr==15)
           a1=axestext_c(x_lab,y_lab, [' '] );  
     else  
           a1=axestext_c(x_lab,y_lab, ['(',label_1,')'] );
     end
     
        if (iso_nr==9 || iso_nr==14)
         label_size=26;
        else
         label_size=20;
        end    

    set(a1,'FontWeight','bold','FontSize',label_size,'rotation',90);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    if show_colorbar==1
    set(br, 'FontSize',18, 'FontWeight', 'bold'); 
    end    
    %%%%%%%%%%%%
    % Grid
     gridm('GLineStyle','--','Gcolor',[.1 .1 .1],'GLineWidth',1.5,...
    'MLineLimit',[lat1 lat2],'Galtitude', .02)

%%%%%%%%%%%%%%%%%%%%%%
% Eastern Ross Sea box for SIC trend
eRS_box=1;         
if iso_nr==9 && param_nr==2 && eRS_box==1 
        
        lon_b_c=[-140 -140 -166 -166 -140];
        lat_b_c=[ -74 -70 -70 -74 -74];
        plotm(lat_b_c,lon_b_c,'-','LineWidth',3,'Color',[0.9 .1 .9])           
end

%%%%% save figure %%%%%%%%
if iso_nr==7 ||iso_nr==11 || iso_nr==12  || iso_nr==14 || iso_nr==15 
    filename=['HadISST_regression_',name,'_',pc_str,'_',season,'_',type_str,'_',num2str(yr_s),'-',num2str(yr_e)];

elseif iso_nr==9
        iso='trend';
    filename=['HadISST_regression_',name,'_',pc_str,'_',season,'_',type_str,'_',num2str(yr_s),'-',num2str(yr_e)];
end

filedir ='C:\PHD\matlab_storage_of_output_files\figures\';
savefilename_c=strcat(filedir,filename);

 % save as png 
orient landscape
% increase quality of picture by increase to -r500

quality_level_nr=3;%%%%%%%%%%%%%%%%%%%%%%%%%

if quality_level_nr==1
   quality_level_str= '-r110';
elseif quality_level_nr==2
  % quality_level_str= '-r150';
   quality_level_str= '-r190';
elseif quality_level_nr==3   
   quality_level_str= '-r240';
end

  export_fig('-png','-nocrop','-painters', '-depsc','-opengl', quality_level_str, savefilename_c); % PNG-nocrop'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% for saving PSA pattern surfaces (step 1)-= one more step=- 
% save regression surface  
folder_c='C:\PHD\matlab_storage_of_output_files\';

    savefilename =[folder_c,filename,'_',num2str(p_level_rc),'_c.mat'];
 
 save(savefilename,'s_all_c','s_all'); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% IPO SAM version (use)
% Where do SAM and PSA1 overlap, SST and SIC
% PSA2
p_level_rc=0.05;
%p_level_rc=0.01;
 type_nr=2; % (2) monthly
 
 yr_s=1979;yr_e=2014;
 
 sign_nr=1; % (1) SAM*(-1)+PSA1, (2) SAM +PSA1
 
 if sign_nr==1
     fact_nr=-1;
 elseif sign_nr==2
     fact_nr=1;
 end
 
 
 sea_nr=1;
 
 folder_c='C:\PHD\matlab_storage_of_output_files\';
%%%%%%%%%%%%%%%

 season='annual';

%%%%%%%%%%%%%%
if type_nr==1
type_str=[];

elseif type_nr==2
type_str='monthly_';

end

param_nr=2;% (1-SST/ 2 SIC)


if param_nr==1
param_str='SST';    
elseif param_nr==2
    
alt_pc=2;%(1) sam psa1 (2) pcs SIC %%%%%%%    
param_str='SIC';

end

%%%%%%%%%%%%%%%%%%%
      season=[];
      
      if param_nr==1
        load([folder_c,'HadISST_regression_',param_str,'_SAM','_', season,'_',type_str,num2str(p_level_rc),'.mat']);
      elseif param_nr==2
          if alt_pc==1
            load([folder_c,'HadISST_regression_',param_str,'_SAM','_', season,...
                '_',type_str,num2str(yr_s),'-',num2str(yr_e),'_',num2str(p_level_rc),'.mat']);
          elseif alt_pc==2
            load([folder_c,'HadISST_regression_',param_str,'_PC1','_', season,...
                '_',type_str,num2str(yr_s),'-',num2str(yr_e),'_',num2str(p_level_rc),'_c.mat']);  
          end
       end
      
% end
%     

s_all_c_sam=s_all_c;  % this is the p-level surface as rs with lower p-levels are masked above

%%%%%%%%%%%%
 ind_sam_nonnan=~isnan(s_all_c_sam);
%%%%%%%%%%

s_all_c_sam(find(s_all_c_sam<0))=-1;  % Negative values -1
s_all_c_sam(find(s_all_c_sam>0))=1;  %
s_all_c_sam(find(s_all_c_sam==0))=NaN; % % remove zeros outside of the sea ice zone

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PSA1

              season=[];
              
              if param_nr==1
                load([folder_c,'HadISST_regression_SST_PSA1','_', season,'_',type_str,num2str(p_level_rc),'.mat']); 
              elseif param_nr==2
                if alt_pc==1
                    load([folder_c,'HadISST_regression_',param_str,'_PSA1','_', season,'_',type_str,num2str(yr_s),'-',num2str(yr_e),'_',num2str(p_level_rc),'.mat']);
                elseif alt_pc==2
                    load([folder_c,'HadISST_regression_',param_str,'_PC2','_', season,'_',type_str,num2str(yr_s),'-',num2str(yr_e),'_',num2str(p_level_rc),'_c.mat']);     
                end
              end

s_all_c_psa1=s_all_c;

%%%%%%%%%%%%%%%%%
 ind_psa1_nonnan=~isnan(s_all_c_psa1);
%%%%%%%%%%%%%%%%%

s_all_c_psa1(find(s_all_c_psa1<0))=-1;  % Negative values -1
s_all_c_psa1(find(s_all_c_psa1>0))=1;  % Negative values -1
s_all_c_psa1(find(s_all_c_psa1==0))=NaN;

% [~,ind_test_c]=ismember(s_all_c_psa1,s_all_c_psa2); % is not member of 

%%%%%%%%%%%%%%
 ind_psa_sum_nonactive=ind_sam_nonnan+ind_psa1_nonnan;
%%%%%%%%%%%%%
    IPO_nr=2; % (1) IPO+, (2) IPO -  % just a name for SIC
    
    if IPO_nr==1
        fac_i=-1;
        ipo_str='ipo_p_';
    elseif IPO_nr==2
        fac_i=1;
        ipo_str='ipo_n_';
    end

ind_psa_sum=s_all_c_sam*( fac_i)+s_all_c_psa1; % multiply by -1 negative corr
ind_psa_sum_c=ind_psa_sum; 
ind_nan= isnan(ind_psa_sum_c);
ind_psa_sum_c(ind_nan)=9999; 
%%%%%%%%%%%%%%%%%%% 
%% save
if sea_nr==1
    if param_nr==1
        filename=['ind_sam_psa1_sum_',ipo_str,num2str(p_level_rc),'.mat'];
    elseif param_nr==2
        if alt_pc==1    
            filename=['ind_',param_str,'sam_psa1_sum_',ipo_str,num2str(p_level_rc),'.mat'];
        elseif alt_pc==2 
            filename=['ind_',param_str,'pc1_pc2_sum_',ipo_str,num2str(p_level_rc),'.mat'];
        end
    end

elseif sea_nr>=2
filename=['ind_sam_psa1_sum_',season,'_',num2str(p_level_rc),'_c.mat'];
end

 savefilename =[folder_c,filename]; 
% save(savefilename,'ind_psa_sum'); 
 save(savefilename,'ind_psa_sum_c','ind_psa_sum_nonactive'); 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%
% Regression ERA-Interim
%    
% Regression_spatial_ERA_I_z500_SAT_DE18_c2.m
% Daniel Emanuelsson 2018
% Matlab 2017a
% Github version 1
%%%%%%%%%%
% subfunctions
% * fig.m fileexchange
% * export_fig.m fileexchange
%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

site='RICE';

%%%%%%%%%%%%%%%%%%%%%%%%%%%
cor_nr=4;  % (1) just trend (4) PCs from EOF z500, (6) PCs from EOF 2mT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

type_nr=2; % (1) annual (2) monthly [use]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
era_name_nr=5; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% 1. z500 **** (Fig 2)
% 3. v850 **** (Fig. 11c) trend
% 5. 2mT SAT ***** (Fig. 3)

if era_name_nr==1
    name='z500';
elseif era_name_nr==3 % uesd in Fig. 11c
    name='v850';
elseif era_name_nr==5
    name='2mT'; 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        yr_nr_s=1;%%%%%%%%%%%
        
        if yr_nr_s==1
            yr_s=1979;
            yr_e=2014;
        elseif yr_nr_s==2 % uesd in Fig. 11c
            yr_s=2000;
            yr_e=2013;    
        end
        
    

        sea_nr=1; %%%%%% (1) Annual 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
%         if sea_nr==1
            season='annual';
%         end


%%%%%%%%%%Figure format
 figure_format=2; % (1) - EPS, (2) - PNG 
%%%%%%%%%%
proj='stereo';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read in ERA-interim data
 addpath C:\PHD\ERA_interim
 
%   ncdisp('ERA_int_monthly_SST.nc')

name_c='ERA_int_monthly_z500_2.nc';
%   ncdisp(name_c)

era_time=ncread(name_c,'time');
era_long=ncread(name_c,'longitude');
era_lat=ncread(name_c,'latitude');
% era_sst=ncread(name_c,'sst');
clear name_c

% era_date=((double(era_time)./24)./365.2400)+1900;

dayssince111=era_time/24;
datevalue=dayssince111+datenum(1900,1,1);
 date_vec=datevec(double(datevalue)); 
 yyyy=date_vec(:,1); 
 mm=date_vec(:,2);
 era_year_num=yyyy+(mm/12-1/12); 
 
%%%%%%%%%%%%%%%%%%%%%%%
 

if strcmp(name,'2mT')==1
    % ECMWF ERA-Inerim (Dee et al. 2011)
    % Monthly means of daily means
    % surface 2mT
    
    name_c='ERA_int_monthly_2m_T.nc';
    %   ncdisp(name_c)
    era_T=ncread(name_c,'t2m'); % 2m temp
    era_T=era_T- 273.15;
    letter='a';
    lat1=-90;
    lat2=-20; 
    label_1='°C s.d.^-^1';
 


 
 elseif strcmp(name,'z500')==1 && strcmp(season,'annual')==1
    % ECMWF ERA-Inerim (Dee et al. 2011)
    % Monthly means of daily means
    % geopotentail 500 hPa
    
    name_c='ERA_int_monthly_z500_2.nc'; 
    era_z500=ncread(name_c,'z'); % GPH z500 %%%%%% file _2 with 2014 and 2015 too ......
    era_z500=era_z500/9.80665;
    letter='a';
    label_1='m';
  
        if strcmp( proj,'stereo')==1
            lat1=-90;
            lat2=0; 
  
        elseif strcmp( proj,'mercator')==1
            lat1=-90;
            lat2=40;
            lon1=-260;
            lon2=-30;
   
        end

  
 elseif strcmp(name,'v850')==1
    % ECMWF ERA-Inerim (Dee et al. 2011)
    % Monthly means of dail means
    % V wind, pressure level 850 hPa 
    
    name_c='ERA_int_monthly_v850.nc'; 
    era_v850=ncread(name_c,'v') ; % v850 'm s**-1'
    letter='c';
    lat1=-90; 
    lat2=-60;  
    if type_nr==1
    label_1='m s^-^1';
    end
 
 end


era_count=find(era_year_num==yr_e)+11;  % data in months
%era_start=find(era_year_num==1980-1);

%if strcmp(season,'annual')==1 % annual 
era_start=find(era_year_num==yr_s); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mm_in_c=[era_start:era_count]';

%%
 %ERA_M_annual_reg
 if cor_nr==1 
     
         iso=[]; 
     
        if      type_nr==1 % annual
 
            wr_yr=[1979:2014]';
            yr_1=find(wr_yr==yr_s);       
            yr_2=find(wr_yr==yr_e);
            
                if sea_nr==1 % annual
                y=wr_yr(yr_1:yr_2);
                elseif sea_nr>=2  % Seasonal (1980-
                y=wr_yr(yr_1+1:yr_2);
                end
                
        elseif   type_nr==2 % monthly
                
                y=[1:era_count]';
            
%                 if sea_nr==1
                y=y(mm_in_c);
%                 elseif sea_nr>=2

                
        end
                
                
            label_2='yr^-^1';
                if era_name_nr==1
                label_1='m decade^-^1';
                elseif era_name_nr==3
                label_1='m s^-^1 decade^-^1';
                end

         
                 letter='c';
          
  elseif cor_nr==4        
         

              load('C:\PHD\matlab_storage_of_output_files\ERA-Interim_PCs_z500_lim0-360_-20_-90_1979-2014_annual');           
            
                  label_1='m s.d.^-^1';
    
    
        if type_nr==1   %%%%%%%%% Annual 
            time_c=MA_PCs_save(:,1);
            
            
        elseif type_nr==2 %%%%%%%%%%% Monthly
            time_c=MA_PCs_save_monthly(:,1);
        end
        
        
        
            yr_1=find( time_c==yr_s);       
            yr_2=find( time_c==yr_e);
    
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            PC_nr=2; % (2) SAM (3) PSA1 (4) PSA2
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% (1) is just time
            
            
            
            label_2='';
            iso='z500';
                                      paper_nr=2; % (1) iso (2) EOF 
    
        if PC_nr==2  % SAM
            
                fact_nr=1;  %2014
      
                 if sea_nr==1
                                     
                    if era_name_nr==1
                            variance_exp=23;
                            letter='a';
                    elseif era_name_nr==3
                            letter='c';
                    elseif era_name_nr==5
                            letter='d';
                    end  
                end
                
            
 
            elseif PC_nr==3 % PSA1
                fact_nr=1; %2014
                if sea_nr==1
                    variance_exp=11;
                        if era_name_nr==1
                            letter='b';
                         elseif era_name_nr==5
                            letter='e';
                           % letter='f';
                         end
                end
                
                
    
              
            elseif PC_nr==4 % PSA2
                 fact_nr=-1;  %2014  To get sign as in Kidson 1988 paper
                 if sea_nr==1
                     variance_exp=10;
                  
                        if era_name_nr==1
                            letter='c';
                        elseif era_name_nr==5
                             letter='f';
%                            letter='e';

                        end 
                 end
        end
            
            
            
         if type_nr==2 %  monthly % works for all seasons
%                 if sea_nr==1
               % y=MA_PCs_save_monthly((yr_1:yr_2),PC_nr)*fact;%%%%%%%%%%%%%%%%%%%%%%%%%%%
               
              y=MA_PCs_save_monthly((mm_in_c),PC_nr)*fact_nr;
                
            %    y=MA_PCs_save_monthly(:,PC_nr)*fact_nr;
                                 
         end
  
 
        
 %%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 elseif cor_nr==6   % SAT PCs      

       load('C:\PHD\matlab_storage_of_output_files\ERA-Interim_PCs_2mT_lim0-360_-30_-90_1979-2014_annual.mat');
       
             if era_name_nr==6
                  label_1='C s.d.^-^1';
             end
                   
    
    
        if type_nr==1   %%%%%%%%% Annual 
            time_c=MA_PCs_save(:,1);
            
            
        elseif type_nr==2 %%%%%%%%%%% Monthly
            time_c=MA_PCs_save_monthly(:,1);
        end
           
        
            yr_1=find( time_c==yr_s);       
            yr_2=find( time_c==yr_e);
    
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            PC_nr=4; % (2) EOF1 (3) EOF2 (4) EOF3
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
       
                label_2='';
                iso='2mT';
                
             if PC_nr==2  % EOF1
                  fact_nr=1;  
                  
                 if sea_nr==1
                 
                    if era_name_nr==1    
                        letter='b';
                    elseif era_name_nr==2 &&  cor_nr==4
                        letter='b';
                    elseif era_name_nr==2 &&  cor_nr==1
                        letter='a';
                    elseif era_name_nr==5
                        letter='a';
                        variance_exp=20; 
                    end
                 
                 end
                
            
 
            elseif PC_nr==3
                   fact_nr=1;  % EOF2
                   
                if sea_nr==1
                    
                        if era_name_nr==1
                            letter='c';
                        elseif era_name_nr==5 
                            letter='c';
                            variance_exp=13;
                        end 
                end
                
    
              
            elseif PC_nr==4 %EOF3
                 fact_nr=1;  
                 
                 if sea_nr==1
                  
                        if era_name_nr==1
                            letter='d';
%                         elseif era_name_nr==2
%                             letter='d';
                        elseif era_name_nr==5
                            letter='b';
                            variance_exp=8;
                        end         
                 end
            end
             
                
            %    y=MA_PCs_save_monthly(:,PC_nr)*fact_nr;
              y=MA_PCs_save_monthly((mm_in_c),PC_nr)*fact_nr;  
                                 
 end
 
%% 
%%%%%%%%%%%%%%%%%%%
% regression of deseasonal 
%%%%%%%%%%%%%%%%%

% method_nr=2;  % Use (2) regstats. It works fine

  
    
if type_nr==2  && era_name_nr==1  %%%%% Monthly data
    
    load('era_z500_c4.mat'); % cosweight and monthly clim removed global z500 not downscaled (run start of eof code)
    %eof_code_run_c4_annual_seasonal.m
    
    ntim=444;
    nlat=241;
    nlon=480;
    
    era_z500_c6_reg_c2 = reshape(era_z500_c4,ntim, nlat,nlon); 
    era_z500_c7_reg_c2=permute(era_z500_c6_reg_c2,[3 2 1]);
    
    %X=era_z500_c7_reg_c2(:,:,(yr_1:yr_2));%%%%%%%%%%%%%%%%%%%%%%%%
     X=era_z500_c7_reg_c2(:,:,(mm_in_c));
     
  %%%%%%%%%%      v850
     elseif type_nr==2  && era_name_nr==3%%%%% v850 Monthly data  
         
     filedir ='C:\PHD\matlab_storage_of_output_files\';
     
      load([filedir,'ERA_interim_v850_clim.mat']); % cosweight and monthly clim removed global z500 not downscaled (run start of eof code)
    %eof_code_run_c4_annual_seasonal.m
    
%         ntim=444; % already done for this one
%         nlat=241;
%         nlon=480;
%     
%         era_z500_c6_reg_c2 = reshape(era_z500_c4,ntim, nlat,nlon); 
%         era_z500_c7_reg_c2=permute(era_z500_c6_reg_c2,[3 2 1]);
    
    %X=era_z500_c7_reg_c2(:,:,(yr_1:yr_2));%%%%%%%%%%%%%%%%%%%%%%%%
        X=ERA_interim_v850_clim(:,:,(mm_in_c));    
        
        
  %%%%%%%%%%      
   
  elseif type_nr==2  && era_name_nr==5%%%%% 2mT Monthly data    
            filedir ='C:\PHD\matlab_storage_of_output_files\';
            load([filedir,'era_2mT_c8.mat']); % cosweight and monthly clim removed global z500 not downscaled (run start of eof code)     
            
        ntim=444;
        nlat=241;
        nlon=480;
            
            era_z500_c6_reg_c2= reshape(era_z500_c8,ntim, nlat,nlon); 
            era_z500_c7_reg_c2= permute(era_z500_c6_reg_c2,[3 2 1]);
    
    %X=era_z500_c7_reg_c2(:,:,(yr_1:yr_2));%%%%%%%%%%%%%%%%%%%%%%%%
            X=era_z500_c7_reg_c2(:,:,(mm_in_c));      
     
    
end
    
[A B C]=size(X);
% b=zeros(2, A,B);bINT=zeros(2,2,A,B);r=zeros(324,A,B);rINT=zeros(324,2,A,B);
b=zeros(2, A,B);
bINT=zeros(2,2,A,B);
r=zeros(C,A,B); % residual
rINT=zeros(C,2,A,B);
STATS=zeros(1,4,A,B);  % 3 p-value

 p_all=zeros(A,B,1);  s_all=zeros(A,B,1);  % just one layer needed?


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% Regression regstats function
 
for i=1:A
    for j=1:B
                        
            s=regstats(squeeze(X(i,j,:)),y,'linear','all');
            s.tstat;
            s_all(i,j,:)=s.tstat.beta(2); % beta- regresion coeficient 
            p_iso_rec=s.tstat.pval;
            p_all(i,j,:)= p_iso_rec(2); % The second element is the p-value for a significantly non-zero slope.
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
%          [a, MSGID] = lastwarn(); % Comment/uncomment if needed
%          warning('off', MSGID)     
            
    end
end

% turn off warning
% w = warning('query','last')
% id = w.identifier
% warning('off',id)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
corr_label=1; %(1/0) %%%%%%%%%%% show colorbar and label
lock_scalebar=1; % (1/0)
%%%%


if type_nr==2 && cor_nr==1%%% TREND only    Monthly
    s_all_c=s_all*10*12; % per decade. %%%%%%  Agree with Bromwich et al. 2012
    

elseif   type_nr==2 || type_nr==1
   s_all_c=s_all;  % unit per s.d.   
end
%%%%%%%%%%%%%%%%%%%%%%%

[c_max, c_min, c_limit]=Ma_max_min(s_all_c);
c_max_c=c_max;
c_limit_c=c_limit;

    if lock_scalebar==1
             
         if    era_name_nr==5 && (cor_nr==4 || cor_nr==6)
%                  c_limit= 1.7879;% 2mT EOF1 2014 max                  
                   c_limit= 1.51;% 2mT EOF1 2014
                                
         elseif    era_name_nr==1 && cor_nr==4 && sea_nr==1
                         c_limit= 47.5260;
                         
         elseif    era_name_nr==3 && cor_nr==1 && yr_nr_s==2         
                         
                        c_limit= 1.00; % 1.6590
         end
      
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fig

fig('units','inches','width',10,'height',9,'font','Helvetica','fontsize',16,'border','on');

   axesm( proj,'MapLatLimit',[lat1 lat2],'Grid','on','ParallelLabel','on','Frame','on',... %
       'MeridianLabel','on','FontWeight','bold','FontSize',18,...  % font size of degree labels
       'mlabellocation',[0:30:180,0:-30:-180]); 

  
set(gca,'box','off','XColor',[1,1,1],'YColor',[1,1,1]);  % turns off border box and axes  
hSurf=surfm(double(era_lat),double(era_long),squeeze(s_all_c(:,:,1))');
hold on
 
        if corr_label==1
            br=colorbar('FontWeight', 'bold', 'FontSize',18);   
            pos=get(br,'Position');  
        
            if cor_nr==4
            
            pos(1)=pos(1)+0.06;      %%%%%%%%%%%%%%%%%%%% x-position of colorbar 
            pos(2)=pos(2)+0.15;      
            pos(3)=pos(3)- 0.0045;  % widthcolorbar
            pos(4)=pos(4)-0.3;  % height colorbar 

            elseif cor_nr==1 && era_name_nr==3
            
            pos(1)=pos(1)+0.02;      %%%%%%%%%%%%%%%%%%%% x-position of colorbar 
            pos(2)=pos(2)+0.09;  
            pos(3)=pos(3)+ 0.000;  % widthcolorbar
            pos(4)=pos(4)-0.195;  % height colorbar  
                    
            else
            pos(1)=pos(1)+0.06;      %%%%%%%%%%%%%%%%%%%% x-position of colorbar 
            pos(2)=pos(2)+0.15;      
            pos(3)=pos(3)- 0.0045;  % widthcolorbar
            pos(4)=pos(4)-0.3;  % height colorbar    
        
            end
        
        set(br,'Position',pos)     
        
        end
        
    
      colormap(b2r(-c_limit,c_limit));
      
      
      if era_name_nr==1 % z500 %%%%%%%%
        color_alt=4;
      elseif era_name_nr==5
        color_alt=1;
      else
        color_alt=1;    
      end
      
      
        if color_alt==1
            colormap(brewermap(256,'*RdBu'));
        elseif color_alt==2
            colormap(flipud(cbrewer('div','Spectral',10)));
        elseif color_alt==3
            colormap(flipud(cbrewer('div','PiYG',20)));              
        elseif color_alt==4
            colormap(flipud(cbrewer('div','RdYlGn',20)));          
        elseif color_alt==5    
             colormap(cbrewer('div','BrBG',20));
        end
      
      
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % p=0.05 significance level contour      
      
    if cor_nr==4 || cor_nr==5   || cor_nr==6 || cor_nr==3 || cor_nr==1    
        co_pa=[.1 .1 .1];
        
    else
        co_pa=[.9 .9 .9];
    end


    h1= contourm( double(era_lat),double(era_long),squeeze(p_all(:,:,1))', [0.05],'--','ShowText','off',...
    'Linecolor',co_pa,'LineWidth', 2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
coast_nr=0;%%%%(1/0) on/off %%%%% Coast lines %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if coast_nr==1

    color_nr=1;
    if color_nr==1
       color_code=[.1 .1 .1];   
    elseif color_nr==2
       color_code=[.6 .6 .6];  
    end

load coast
%%%%%%%%%%%%%
% to be able to use coastline for continents in combination with ant
% coastline and grounding line from bedmap
in_c=find(lat<-60);
lat_cr=lat;
lon_cr=long;
lat_cr(in_c)=NaN;
lon_cr(in_c)=NaN;

%%%%%%%%%%%%%
plot3m(lat_cr,lon_cr,'-k','LineWidth', 1.5,'color',color_code);
%plot3m(lat(in_c),long(in_c),'.k','MarkerSize', 1);

addpath C:\PHD\matlab_mapping
bedmap2('gl','k','LineWidth', 1.5,'color',color_code);
bedmap2('patchshelves','LineWidth', 1.0,'facecolor','none','frame','on','edgecolor',color_code)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% label above colorbar
label_size=22;
  
if cor_nr==1 && era_name_nr==3 %
    a1=axestext_c(1.101, +0.33, ['(',label_1,')'] );    
elseif cor_nr==4 || cor_nr==5 || cor_nr==1 || cor_nr==6
    a1=axestext_c(0.04, +0.028, ['(',label_1,')'] );    
else
    a1=axestext_c(0.1, -0.04, ['(',label_1,' ', label_2,')'] );
end

set(a1,'FontWeight','bold','FontSize',label_size,'rotation',90);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variance explained by EOF
        if ((cor_nr==4 && era_name_nr==1) || (cor_nr==6 && era_name_nr==5)) && sea_nr==1  
            a2=axestext_c(0.95, +0.98, [num2str(variance_exp),'%'] );
            set(a2,'FontWeight','bold','FontSize',label_size+4);           
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lat long labels. Long  at this lat and placement of lat         
    if cor_nr==6 || cor_nr==4 && era_name_nr==5
        mlab_par=-28;
    elseif cor_nr==1 && era_name_nr==3    
        mlab_par=-64; 
    else
        mlab_par=-10;
    end
        
mlabel('MLabelParallel',mlab_par,'PLabelLocation',[-75 -60 -45 -30 -15  0],'PLabelMeridian',100) 


%%%%%%%%%%%%%%%%%%%%%%
%  Title

if cor_nr==1    
        if type_nr==2
            type_str='Monthly';
        else
            type_str='Annual';    
        end
        
        t4=title(['Trend ', name,' ', type_str,' ',num2str(yr_s),'-',num2str(yr_e)]);   
            
elseif cor_nr==3
        t4= title(['Regression ', iso,'/', name]);
   
elseif cor_nr==4 || cor_nr==5  % EOFs 
    
        if type_nr==2
            type_str='Monthly';
        else
            type_str='Annual';    
        end
    
    if PC_nr==2
        pc_str='SAM';
    elseif PC_nr==3
        pc_str='PSA1';
     elseif PC_nr==4
        pc_str='PSA2';        
    end 
    
    if sea_nr==1
        season='Annual';
    end
     
            if cor_nr==4 && sea_nr==1 &&  era_name_nr==1
                t4=title([ pc_str,' ','EOF',num2str(PC_nr-1)]);   
%             elseif cor_nr==4 && sea_nr==1 &&  era_name_nr==2
%                 t4=title([ pc_str,' ','EOF',num2str(PC_nr-1),' ',type_str]);   
             
             elseif cor_nr==4 && sea_nr==1 &&   era_name_nr==5
                 
%                     if PC_nr<=3
                        t4=title([ pc_str,' ','PC',num2str(PC_nr-1),', SAT ',type_str]);   
%                     elseif PC_nr==4
%                         t4=title([ pc_str,' ','PC',num2str(PC_nr-1),', SAT ',type_str,' (-1)']);       
%                     end
                        
             
            elseif cor_nr==4 && sea_nr>=1 &&  ( era_name_nr==1 || era_name_nr==2 )
            t4=title([ pc_str,' ',name,' ',season]);  
            
%             elseif cor_nr==4 && sea_nr==1 &&  era_name_nr==2
%               t4=title([ pc_str,' ',season,' ', name]);  
                
            elseif cor_nr==5
             t4=title([ pc_str]);      
            end
            
 elseif cor_nr==6           
  
        if type_nr==2
        type_str='Monthly';
        else
        type_str='Annual';    
        end
        
            if PC_nr==2
                pc_str='PC1';
                %pc_str='EOF1';
                
            elseif PC_nr==3
              pc_str='PC2';
              %  pc_str='EOF2';
            elseif PC_nr==4
                pc_str='PC3';
                 %    pc_str='EOF3';
            end

                    t4=title(['EOF',num2str(PC_nr-1),', SAT ',type_str]);  
       
 
end


        pos=get(t4,'Position');
               
% move title        
    if cor_nr==1
        pos(2)=pos(2)-0.030;
    elseif cor_nr==4   
        pos(2)=pos(2)-0.055;
    else    
        pos(2)=pos(2)-0.055;
    end
        set(t4,'Position',pos,'FontSize',22, 'FontWeight', 'bold');  


   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Letter
        letter_pos= [170 710 50 62];
        letter_size=38;
    

          TextBox = uicontrol('style','text');
          set(TextBox,'String',letter,'position',letter_pos,'FontWeight','bold','FontSize',letter_size ); % x position y position xsize ysize
          set(TextBox,'foregroundcolor', [0 0 0], ...
         'backgroundcolor', [1 1 1]);
     
     %%%%%%%%%%%%%%%%%
     % colorbar fontsize 
     if cor_nr==4
        scalebar_fontsize=22;    
     else
        scalebar_fontsize=16;
     end
     
      set(br, 'FontSize',scalebar_fontsize, 'FontWeight', 'bold'); 
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
     % save fig
        str_end='_c';
     
     
        if cor_nr==4 || cor_nr==6      
            filename=['ERA_interim_regression_',name,'_',pc_str,'_',season,...
                '_',type_str,'_',num2str(yr_s),'-',num2str(yr_e),str_end];
        
        elseif cor_nr==1 && era_name_nr==3
            filename=['ERA_interim_regression_',name,'_',season,'_',type_str,'_',num2str(yr_s),'-',num2str(yr_e)];
        
        else
            filename=['ERA_interim_regression_',name,'_',iso,'_',season];
        end
        

filedir ='C:\PHD\matlab_storage_of_output_files\figures\';

if corr_label==1
savefilename_c=strcat(filedir,filename,'_crop');
else
savefilename_c=strcat(filedir,filename);

end


% save as png 
orient landscape
export_fig('-png','-nocrop','-painters', '-depsc','-opengl', '-r240',savefilename_c); % PNG-nocrop'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save EOF surfaces (use the next one)

filedir ='C:\PHD\matlab_storage_of_output_files\';

      filename_c=['ERA_interim_regression_surface_',name,'_',pc_str,'_',season,'_',type_str,'_',num2str(yr_s),'_',num2str(yr_e)];
      savefilename_c=strcat(filedir,filename_c);

      %save

    save(savefilename_c,'s_all','p_all'); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%% For saving regression pattern surfaces -= one more step=- (use)
 %  s_all_c=s_all_c1;

p_level_rc=0.05;
%p_level_rc=0.01;

s_all_c(p_all>=p_level_rc)=NaN; 


%%%%%%%%%%%%%
% [wa, wb]=size(s_all_c); 
% save regression surface  
folder_c='C:\PHD\matlab_storage_of_output_files\';

savefilename =[folder_c,filename,'_',num2str(p_level_rc),'_c.mat'];

save(savefilename,'s_all_c'); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
 %% Alt 2 save surface -= Loads saved surfaces from above =-
 
p_level_rc=0.05;

%  
  folder_c='C:\PHD\matlab_storage_of_output_files\';

    if era_name_nr==5
        load([folder_c,'ERA_interim_regression_2mT_PC1_annual_Monthly_',num2str(p_level_rc),'_c.mat']);
    elseif era_name_nr==1
        load([folder_c,'ERA_interim_regression_z500_SAM_Annual_Monthly_',...
            num2str(yr_s),'-',num2str(yr_e),'_c_',num2str(p_level_rc),'_c.mat']);

    end      

    

s_all_c_pc1=s_all_c;  % this is the p-level surface as rs with lower p-levels are masked above

%%%%%%%%%%%%
 ind_pc1_nonnan=~isnan(s_all_c_pc1);
%%%%%%%%%%

s_all_c_pc1(find(s_all_c_pc1<0))=-1; % Negative values -1
s_all_c_pc1(find(s_all_c_pc1>0))=1;  % positive values 1


    if era_name_nr==5
            load([folder_c,'ERA_interim_regression_2mT_PC3_annual_Monthly_',num2str(p_level_rc),'_c.mat']);
    elseif era_name_nr==1
            load([folder_c,'ERA_interim_regression_z500_PSA1_Annual_Monthly_',...
            num2str(yr_s),'-',num2str(yr_e),'_c_',num2str(p_level_rc),'_c.mat']);
    end

s_all_c_pc3=s_all_c;

%%%%%%%%%%%%%%%%%
 ind_pc3_nonnan=~isnan(s_all_c_pc3);
%%%%%%%%%%%%%%%%%

s_all_c_pc3(find(s_all_c_pc3<0))=-1;  % Negative values -1
s_all_c_pc3(find(s_all_c_pc3>0))=1;  % positive values 1



% [~,ind_test_c]=ismember(s_all_c_psa1,s_all_c_psa2); % is not member of 

%%%%%%%%%%%%%%
 ind_psa_sum_nonactive=ind_pc1_nonnan+ind_pc3_nonnan;
%%%%%%%%%%%%%

 ind_psa_sum=s_all_c_pc1+s_all_c_pc3;
 ind_psa_sum_c=ind_psa_sum;
 
ind_nan= isnan(ind_psa_sum_c); 
ind_psa_sum_c(ind_nan)=9999; 
 


    if era_name_nr==5
        filename=['pc_2mT_ind_psa_sum_',num2str(p_level_rc),'.mat'];
    elseif era_name_nr==1
        filename=['pc_z500_ind_psa_sum_',num2str(p_level_rc),'.mat'];
    end


 savefilename =[folder_c,filename]; 
% save(savefilename,'ind_psa_sum'); 
 save(savefilename,'ind_psa_sum_c','ind_psa_sum_nonactive'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

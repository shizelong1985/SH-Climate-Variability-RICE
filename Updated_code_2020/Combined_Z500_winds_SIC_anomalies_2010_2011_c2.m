%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% combined_fields_figure_Palmer_corr_SIC_z500_winds.m
% Combine SIC, Z500, Winds anoamly fields for 2010-2011
% Daniel Emanuelsson, 2020
% MATLAB 2018a 
% Github version 1 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
clear all
close all
    step_nr=5;
    name='SIC';
    % saved from Anomalies_HadISST_SIC_c2.m
load(['C:\PHD\matlab_storage_of_output_files\anomalies_2010_2011_SIC.mat'])

% Defaults
  proj='stereo';
  lat1=-90;
  lat2=-50;
%  lat2=-40; 
%  lat2=-10;  
  box_use=0;
  letter='b';
 
  
  lock_scalebar=1;%%%%%%%%%
  season_str='Annual';
 %%%%%%%

  f2=fig('units','inches','width',10.5,'height',10.5,'font','Helvetica','fontsize',16,'border','on'); 

  axesm( 'stereo','MapLatLimit',[lat1 lat2],'Grid','on',...
      'Frame','on','ParallelLabel','on',...
       'MeridianLabel','on','FontWeight','bold','FontSize',18, 'mlabellocation',[0:-30:-360, 0:30:180]); 
 
  set(gca,'box','off','XColor',[1,1,1],'YColor',[1,1,1]);  % turns off border box and axes 

  hold on

  set(gca,...
     'linewidth',2,...
     'FontWeight','bold' );
      
 gridm('GLineStyle','--','Gcolor',[.5 .5 .5],'GLineWidth',1.5,... %%%%%%%%%% Grid
    'MLineLimit',[lat1 lat2],'Galtitude', .02)



  
if lock_scalebar==1 
   
     c_limit=0.08; % 
         
end


% 

% Rc6=squeeze(Rc5(1,2,:,:))';
% Rc6(:,361)=Rc6(:,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if type_nr==1
%     wr_c=s_all*10;
% elseif  type_nr==2   
%     wr_c=s_all*10*12;
% end

hSurf=surfm(double(HadISST_lat),double(HadISST_lon),X_surface');
colormap(b2r(-c_limit,c_limit));

hold on
% colormap(brewermap(256,'*RdBu'))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    color_alt=6;    



if color_alt==1
    colormap(brewermap(256,'*RdBu'));
elseif color_alt==2
    colormap(flipud(cbrewer('div','Spectral',10)));
elseif color_alt==3
    colormap(flipud(cbrewer('div','PiYG',20)));      
elseif color_alt==4
    colormap(flipud(cbrewer('div','RdGy',20)));      
elseif color_alt==5
    colormap(flipud(cbrewer('div','RdBu',20)));    
elseif color_alt==6
    colormap(cbrewer('div','BrBG',20));
elseif color_alt==7
    colormap(flipud(cbrewer('div','RdYlBu',20)));    
end

%%%%%%%%%%%%%%%%%%%%%%
% hSurfAx=(gca);
cRange= caxis;
% contourm function keeps projection defined above
%     h1= contourm( ndsic_lat,ndsic_lon,p_all', [0.05],'-k','ShowText','off',...
%    'LineWidth', 2);

% Then reset the color axis according to the range determined above:
caxis(cRange); 
alpha(1) 
%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%

%c_limit_c
contour_lim=1.3;
% max_contour=1; % (1/0)
polygon_nr=2; % (1/2)
% contour_p_level=0.05;

% if max_contour==1
%     contour_lat_lon=contourm(ndsic_lat,ndsic_lon,squeeze(P(1,2,:,:))',...
%        [contour_p_level contour_p_level],'-k','ShowText','off',...
%     'Linecolor',[.1 .1 .1],'LineWidth', 2);%,'LineColor','none');  %%%%%%%%%%%%%%%%%
% end
%%%%%%%%%%%%%%%%%%%%%%%%%
% ind_lon = find(HadISST_lon<=-93.5 & HadISST_lon>=-104);
% ind_lat =find(HadISST_lat<=-33 & HadISST_lat>=-38);
hSurfAx=(gca);
cRange= caxis;

show_colorbar=1;

if show_colorbar==1
  %h=colorbar('SouthOutside'); 
  h=colorbar('EastOutside');  
%   set(h, 'Position', [.76 .145 .015 .75], 'FontSize',18, 'FontWeight', 'bold');  
end

 a1=axestext_c(1.09,+0.53, ['SIC'] ,'rotation',-90,'FontSize',22, 'FontWeight', 'bold');    %050
%     a1=axestext_c(1.09,+0.6, ['Trend (% dec^-^1)'] ,'rotation',-90,'FontSize',22, 'FontWeight', 'bold');
    
    

    % modifiy colorbar
    pos=get(h,'Position');
    
%         pos(1)=pos(1)+0.024;      %%%%%%%%%%%%%%%%%%%% x-position of colorbar   (+) --->>   +0.02;
        pos(1)=pos(1)+0.022; % 0.035
        pos(2)=pos(2)+0.09;  
        pos(3)=pos(3)+ 0.012;  % widthcolorbar
        pos(4)=pos(4)-0.18;  % height colorbar
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
        set(h,'Position',pos)




   mlabel('MLabelParallel',-46  ,'PLabelLocation', [-75 -60 -45 -30],'PLabelMeridian',100) 


%%%%%%%%%%%%%%%%%%%%%%%%%%
% letter
TextBox = uicontrol('style','text');
       if strcmp(name,'SIC')==1 
       letter_position=  [200 820 60 60];%
       else
       letter_position=   [200 820 60 60];
       end
  
%        if yr_e==2009
            set(TextBox,'String',letter,'position',letter_position,'FontWeight','bold','FontSize',44 ); % x position y position xsize ysize 
%        else
%            set(TextBox,'String',letter,'position',letter_position,'FontWeight','bold','FontSize',40 ); % x position y position xsize ysize
%        end
            
% end
            set(TextBox,'foregroundcolor', [0 0 0], ...
            'backgroundcolor', [1 1 1]);
               
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
  %clear all
%     if z_level_nr==1
%        z_level_str='z200';
%     elseif z_level_nr==2
%        z_level_str='z500';  
%     elseif z_level_nr==3
%        z_level_str='z850'; 
%     end

% from Anomalies_era_i_DE20_c2.m  
   load(['C:\PHD\matlab_storage_of_output_files\anomalies_2010_2011_z500_c2.mat'])   

      
%       if season_nr==6
%        yr_s=1979;
%       else
%        yr_s=1980;   
%       end
%       
%     if corr_w_std==1
%     z_level_str=['std ',z_level_str];
%     end    
%       
%     load(['C:\PHD\matlab_storage_of_output_files\Palmer_',iso_str_c2,'_',z_level_str,'_',season,'_ERA_5_correlation _',...
%         num2str(yr_s),'_2011stereo_crop.mat'])      
%      
%   end
  
%   %contourm function keeps projection defined above
%     Pc=squeeze(P(1,2,:,:))';
% if type_nr==1
     Rc_c1=X_surface;
     Rc_cc1=X_surface;
% elseif type_nr==2
%    Rc_c1=s_all*10*12; 
% end

% if type_nr==1
% Rc_cc1=s_all*10;
% elseif type_nr==2
% Rc_cc1=s_all*10*12;    
% end
Rc_cc1(Rc_c1>0)=NaN; % dashed

% Rc_cc1(s_all>0)=NaN;% line contour
% Rc_cc1(p_all>=0.05)=NaN;
% 
% Pc_c_1=p_all;
% Pc_c_1(s_all>0)=NaN;


%%%%%%
% if type_nr==1
Rc_c2=X_surface;
Rc_cc2=X_surface;
% elseif type_nr==2
% Rc_c2=s_all*10*12;    
% end
% Rc_c2(s_all<0)=NaN;


% if type_nr==1
% Rc_cc2=s_all*10;
% elseif type_nr==2
% Rc_cc2=s_all*10*12;    
% end

Rc_cc2(Rc_c2<0)=NaN;
% Rc_cc2(p_all>=0.05)=NaN;
% 
% Pc_c_2=p_all;
% Pc_c_2(s_all<0)=NaN;


     contourm( double(era_lat),double(era_long),Rc_cc1','LevelStep',step_nr, 'LineColor',[0.0 0.2 0.6],'LineStyle',':','ShowText','off','LineWidth', 2); 
  
%      if nanmin(nanmin(Pc_c_1)) <=0.05
%      contourm( double(era_lat),double(era_long),Rc_cc1','LevelStep',step_nr, 'LineColor',[0.0 0.2 0.6],'LineStyle','-','ShowText','off','LineWidth', 2);
%      end
     
     contourm( double(era_lat),double(era_long),Rc_cc2','LevelStep',step_nr, 'LineColor',[1.0 0.1 0],'LineStyle',':','ShowText','off','LineWidth', 2);
 
% not totaly reliable (reason? sometimes to few point need more than one point to make contour), so comment/uncomment
% not working for annual
%      if nanmin(nanmin(Pc_c_2)) <=0.05 %&& season_nr>1 % && ~(season_nr==6 && corr_w_std==0 && movstd_iso==0)
%       contourm( double(era_lat),double(era_long),Rc_cc2','LevelStep',step_nr, 'LineColor',[1.0 0.1 0],'LineStyle','-','ShowText','off','LineWidth', 2);    
%      end
     
     

     
     
     
% hold on

show_wind_corr=1;

%  % wind corr
%  if show_wind_corr==1
%      
%      
%     
%     
%     if w_level_nr==2
%        w_level_str='v500';  
%     elseif w_level_nr==3
%        w_level_str='v850'; 
%     end
%     
%     
%     if corr_w_std==1
%     w_level_str=['std ',w_level_str];
%     end  
    
    
  
%     if  season_nr==1  
    %load(['C:\PHD\matlab_storage_of_output_files\Palmer_',iso,'_',w_level_str,'_',season,'_ERA_5_correlation _1979_2011stereo_crop.mat'])
    
    load(['C:\PHD\matlab_storage_of_output_files\anomalies_2010_2011_v850_c2.mat'])
    
%     else
      
%       if season_nr==6
%        yr_s=1979;
%       else
%        yr_s=1980;   
%       end
%         
%     load(['C:\PHD\matlab_storage_of_output_files\Palmer_',iso_str_c2,'_',w_level_str,'_',season,'_ERA_5_correlation _',...
%         num2str(yr_s),'_2011stereo_crop.mat'])    
%     
% %     end
%   
%     clear w_level_str
    
%     if type_nr==1
    Rc_cv=X_surface;
%     elseif type_nr==2
%     Rc_cv=s_all*10*12;    
%     end 
%     Pc_v=p_all;
    

%     if w_level_nr==2
%        w_level_str='u500';  
%     elseif w_level_nr==3
%        w_level_str='u850'; 
%     end  
%     
%     if corr_w_std==1
%     w_level_str=['std ',w_level_str];
%     end
%     
%     
     clear  X_surface
  
%      if  season_nr==1
     %load(['C:\PHD\matlab_storage_of_output_files\Palmer_',iso,'_',w_level_str,'_',season,'_ERA_5_correlation _1979_2011stereo_crop.mat']) 
     load(['C:\PHD\matlab_storage_of_output_files\anomalies_2010_2011_u850_c2.mat'])
     
%      else
%      
%       if season_nr==6
%        yr_s=1979;
%       else
%        yr_s=1980;   
%       end
%          
%          
%      load(['C:\PHD\matlab_storage_of_output_files\Palmer_',iso_str_c2,'_',w_level_str,'_',season,'_ERA_5_correlation _',...
%          num2str(yr_s),'_2011stereo_crop.mat'])     
%      
%      end
  
%   if type_nr==1
    Rc_cu=X_surface;
%   elseif type_nr==2
%     Rc_cu=s_all*10*12;  
%   end
%   
  
%    Pc_u=p_all;
    
    %%%%%%%%%% Only show significant wind r
    
    [a b c]=size(Rc_cv);
    
%     ind_w_v=zeros(a,b);
%     ind_w_u=zeros(a,b);
% %     ind_w_tot=NaN(a,b);
%     
%     for i=1:a
%        for j=1:b
%         
%            if Pc_v(i,j)<0.1
%                ind_w_v(i,j)=1;              
%            end
%            
%            if Pc_u(i,j)<0.1
%                ind_w_u(i,j)=1;              
%            end          
%                    
%        end
%     end
    
    
%     ind_w_tot=ind_w_v+ind_w_u;
%     
%      Rc_cu(ind_w_tot==0)=NaN;
%      Rc_cv(ind_w_tot==0)=NaN; 
    %     ind_v=find( Rc_cv(Pc>=0.1)); 
       
    
    %%%%%%%%
    
    
    
    
%     clear  Rc_c
    
    
     [long,latg] = meshgrid(era_long,era_lat );
    
     
     latg_c=latg((141:end),:); % only show winds inside figure limits
     long_c=long((141:end),:);
     Rc_cu_c=Rc_cu';
     Rc_cv_c=Rc_cv';
     
     wind_den=1;
     density_limit=10;
     unitString='ms^-^1';
    
     quivermc_c(double(latg_c), double(long_c), Rc_cu_c((141:end),:), Rc_cv_c((141:end),:),...
    'linewidth',1.5,'reference',wind_den,'density',density_limit,'units',unitString,'color',[.9 .0 .0]); % 'reference','median'
     
%  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 coast_nr=0;  % On/Off
if coast_nr==1
    
    load coast
    %%%%%%
    % to be able to use coastline for continents in combination with ant
    % coastline and grounding line from bedmap
    in_c=find(lat<-60);
    lat_cr=lat;
    lon_cr=long;
    lat_cr(in_c)=NaN;
    lon_cr(in_c)=NaN;

    %%%%%%%%%%%%%
    addpath C:\PHD\matlab_mapping
%         if strcmp(name,'SIC')==1
            color_code=[.4 .4 .4];    
            plot3m(lat_cr,lon_cr,'-','LineWidth', 2,'color',color_code);
            bedmap2('gl','LineWidth', 2,'color',color_code);
%         end

    addpath C:\PHD\matlab_mapping
    bedmap2('gl','LineWidth', 2, 'color',color_code );
    bedmap2('patchshelves','LineWidth', 1.0,'facecolor','none','frame','on','edgecolor',color_code)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

if strcmp(name,'SIC')==1 && show_maximum_point==1     
    plotm(ndsic_lat(min_ind_2,min_ind_1), ndsic_lon(min_ind_2,min_ind_1),'.','color',[.4,.8,.2],'MarkerSize',36);  % max corr point

%          TextBox = uicontrol('style','text');
%          set(TextBox,'String',[' r = ',num2str(round(X_r_max*1000)/1000)],...
%              'position',[310 130 170 50],'FontWeight','bold','FontSize',22 ); % x position y position xsize ysize
%           set(TextBox,'foregroundcolor', [0 0 0], ...
%          'backgroundcolor', [1 1 1]);    
end

%%%%%%%%%%%%%%
% marker for ice-core site
% if  strcmp(site,'RICE')==1 || strcmp(site,'Palmer')==1
%         site_coor=star_coord_WA(site); % RICE marker
%         % Circle Alt
%         plotm(site_coor(1),site_coor(2),'.','MarkerSize',30,'MarkerFaceColor','none','MarkerEdgeColor',[.9,.1,.9]); %,'LineWidth',2.5); % [.4,.6,.1]
%         %%%% Dot Alt
%         %plotm(site_coor(1),site_coor(2),'.','MarkerSize',40,'MarkerFaceColor','none','MarkerEdgeColor',[.9,.1,.9],'LineWidth',2.5);   
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% title
if show_title==1 || show_title==2 || show_title==3 || show_title==4
                                   % SST_dataset_str='ndsic';
                                    SST_dataset_str='NSIDC';   
                                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yr_s=2010;
yr_e=2011;      
            
hp1=  title(['Anomalies SIC, Z500, Winds',...
      ' ', num2str(yr_s),' - ',num2str(yr_e)],'FontSize',20, 'FontWeight', 'bold');
                        
 % move title       
        if strcmp(proj,'stereo')==1  
            pos=get(hp1,'Position');
            pos(2)=pos(2)-0.03;
            set(hp1,'Position',pos,'FontSize',28, 'FontWeight', 'bold');    
   
        end
       
end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if show_colorbar==1
        set(h, 'FontSize',20, 'FontWeight', 'bold'); 
    end
    
     

hold off
%%%%%%%%%%%%%%%%%%%%
% save figure
 filedir ='C:\PHD\matlab_storage_of_output_files\figures\'; 
%     if rcontour_psa==1
%                 filename=strcat(filedir,'ndsic_',name,'_',iso,'_',season,'_',proj,'_',num2str(yr_s),'_',num2str(yr_e),'_wPSA_pattern');
%      if season_nr==1 

                filename=strcat(filedir,'anomalies_SIC_Z500_Winds',...
                   '_',proj,'_',num2str(yr_s),'_',num2str(yr_e));
                
%      elseif season_nr==6           
%                 filename=strcat(filedir,'combined_trend_NSIDC_ERA5_',site,'_',name,'_',season_str,'_',z_level_str,...
%                     '_',w_level_str,'_',proj,'_',num2str(yr_s),'_',num2str(yr_e));                
%                 
%      else
%                 filename=strcat(filedir,'combined_trend_NSIDC_ERA5_',site,'_',name,'_',season_str,'_',z_level_str,...
%                     '_',w_level_str,'_',proj,'_',num2str(yr_s+1),'_',num2str(yr_e));
%      end

savefilename_c=filename;

 % save as png 
orient landscape
% increase quality of picture by increase to -r500
% if figure_format==1
%   export_fig('-eps','-nocrop','-painters', '-depsc','-opengl', '-r100',savefilename_c); % EPS works  func needs ghostscript 
% elseif figure_format==2
  export_fig('-png','-painters', '-depsc','-nocrop','-opengl', '-r190',savefilename_c); % PNG  110
  export_fig('-pdf','-painters', '-depsc','-nocrop','-opengl', '-r190',savefilename_c); % PNG  110 
% end

% end

toc
%%%%%%%%%%
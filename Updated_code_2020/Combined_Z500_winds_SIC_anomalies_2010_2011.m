% combined_fields_figure_Palmer_corr_SIC_z500_winds.m


tic
clear all
close all


     step_nr=5;
%     wind_den=1;
%     end

% if season_nr==1
%     season='annual';
%     
% elseif season_nr==2
%      season='DJF';
%      
% elseif season_nr==3
%      season='MAM';
%      
% elseif season_nr==4
%      season='JJA';
%      
% elseif season_nr==5
%      season='SON';
%      
% elseif season_nr==6
%      season='MAMJJA';    
% end


% w_level_nr=3; % (1) 200 (2) 500 (3) 850       for V U winds
% z_level_nr=1;
% 
%  iso_nr=4;
%   
%     if iso_nr==1
%     iso='dD';
%     elseif iso_nr==2
%     iso='Acc';  
%     elseif iso_nr==3
%     iso='MSA';
%     elseif iso_nr==4
%     iso='d18O';
%     elseif iso_nr==5
%     iso='d-excess';    
%     end

    name='SIC';
%     if corr_w_std==1
%     name=['std ',name];
%     end   
%     
%     if movstd_iso==1
%         iso_str_c=['mstd ',iso]; 
%         iso_str_c2=['mstd_',num2str(movstd_win),'_',iso]; % for file name    
%     else
%         iso_str_c=iso;
%         iso_str_c2=iso;
%     end
    
    
    

% if  season_nr==1
load(['C:\PHD\matlab_storage_of_output_files\anomalies_2010_2011_SIC.mat'])
%load(['C:\PHD\matlab_storage_of_output_files\NSIDC_corr_Palmer_SIC_',iso,'_',season,'_stereo_1979_2011.mat'])
% load(['C:\PHD\matlab_storage_of_output_files\NSIDC_corr_Palmer_SIC_',iso_str_c2,'_',season,'_stereo_1979_2011.mat'])
% if type_nr==1
%    load(['C:\PHD\matlab_storage_of_output_files\NSIDC_Reg_SIC_Annual_annual_means_',...
%        num2str(yr_s),'_',num2str(yr_e),'.mat']) 
% end


% elseif season_nr==6
% load(['C:\PHD\matlab_storage_of_output_files\NSIDC_corr_Palmer_',name,'_',iso_str_c2,'_',season,'_stereo_1979_2011.mat']) 
% else
% %load(['C:\PHD\matlab_storage_of_output_files\ndsic_corr_Palmer_SIC_',iso,'_',season,'_stereo_1980_2011.mat'])
% load(['C:\PHD\matlab_storage_of_output_files\NSIDC_corr_Palmer_',name,'_',iso_str_c2,'_',season,'_stereo_1980_2011.mat']) 
% end



tic
% Defaults
  proj='stereo';
  lat1=-90;
  
  
   lat2=-50;
%  lat2=-40; 
%    lat2=-10; 
  
  
  box_use=0;
  

        letter='b';


%   name='SIC';
%  site='RICE';
 
  show_maximum_point=0;
  show_title=1;
  
  lock_scalebar=1;%%%%%%%%%%%%%%%%%%%%%%%
  
  Int_in=0;
  
%   if Int_in==1
%     RSAS_SIC_box=0;    
%   elseif Int_in==0    
%     RSAS_SIC_box=0; % (1) area defined using RICE corr pattern (2) ADP following Yuan 2004
%   end
  
%   if season_nr==1
     season_str='Annual';
%   elseif season_nr==2
%      season_str='DJF';
%   elseif season_nr==3
%      season_str='MAM';
%   elseif season_nr==4
%      season_str='JJA';
%   elseif season_nr==5
%      season_str='SON';
%   elseif season_nr==6
%      season_str='MAMJJA';    
%   end
  
  
  
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


% Rc=s_all;
% 
% Rc5=s_all;
% % Rc2=R;
% 
% % Rc(P>=p_level_rc)=NaN; 
% % Rc2(P>=0.01)=NaN; % Mask out above or equal R's
% Rc5(p_all>=0.1)=NaN; 


% c_min=squeeze(s_all(1,:,:));
% c_min=c_min(2,:);
% c_min=min(c_min);
% 
% c_max=squeeze(s_all(1,:,:)); 
% c_max=c_max(2,:);
% c_max=max(c_max);
% 
% c_limit=max(abs(c_min),abs(c_max));
% c_limit_c=c_limit;


% if strcmp(name,'SIC')==1
% 
%   X_r=squeeze(Rc(1,2,:,:))';
%   X_r(:,(161:316))=NaN;
%   [X_r_min, min_ind_1]=nanmin(nanmin(X_r)); % negative corr area RS
%   [X_r_min_c, min_ind_2]=nanmin(nanmin(X_r,[],2));
%   
% 
%   Px_c=squeeze(p_all(1,2,:,:))';
%   
%   p_level(Px_c(min_ind_2,min_ind_1)); % p level at correlation extreme

  
if lock_scalebar==1 
   
     c_limit=0.08; % 
         
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%
% if corr_w_std==1
%    X_r=squeeze(Rc(1,2,:,:))';
% 
%   X_r_max=nanmin(X_r,[],2);
%   min_or=[1:332]';
%   
%   Ma_order=NaN(332,1);
%   Ma_order(:,1)=X_r_max;
% %   Ma_order(:,2)=min_or;
%   
%   [Ma_order_c, index_o]=sort(Ma_order,1);
%   %%%%%%%%%%%%
%   % std SIC corr
%   %ord_nr=1;  % noise
%   %ord_nr=2; % Bell S west -0.6096
%   %ord_nr=3; % noise
%   %ord_nr=4; % Bell S west -0.6060, same region as above
%   % ord_nr=5; % noise
%   % ord_nr=6; % noise
%   % ord_nr=7; % noise
%   % ord_nr=8; % Bell S east -0.5964
%   % ord_nr=9; % ||
%    ord_nr=10; % eastern Weddell Sea -0.5898
%   % ord_nr=11; % noise
%   % ord_nr=12; % noise
%    
%   in_or1=index_o(ord_nr);
%   in_or2=find(X_r(index_o(ord_nr),:)==Ma_order_c(ord_nr));
%   
%   in_or3=index_o(2);
%   in_or4=find(X_r(index_o(2),:)==Ma_order_c(2));  
%   
%   % plot min location below
%   
%   % Ma_order_c(ord_nr)
%   
% %   if  ord_nr==2 % Bell
%   color_id1=[1 .2 .2];
% %   elseif ord_nr==10 % Weddell Sea
%   color_id2=[1 .7 .1];
% %   else
% %   color_id=[.4,.8,.2];    
% %   end
% 
% end


%%%%%%%%%%%%%%%% Fill in seam  % only one side
% % HadISST_lat_c=HadISST_lat;
% HadISST_lon_c=HadISST_lon;
% HadISST_lon_c(361,:)=179.99;
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
if Int_in==1 % color for Fig. 10
   color_alt=1;   
else
    color_alt=6;    
end


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

    
% if yr_e==2009 % suppl. mat
 a1=axestext_c(1.09,+0.53, ['SIC'] ,'rotation',-90,'FontSize',22, 'FontWeight', 'bold');    %050
%     a1=axestext_c(1.09,+0.6, ['Trend (% dec^-^1)'] ,'rotation',-90,'FontSize',22, 'FontWeight', 'bold'); 
% else
%     
%     a1=axestext_c(1.065,+0.57, ['Correlation'] ,'rotation',-90,'FontSize',20, 'FontWeight', 'bold');    %050
% 
% end 
    
    



    
    
    % modifiy colorbar
    pos=get(h,'Position');
    
%     if yr_e==2009
    
        
%         pos(1)=pos(1)+0.024;      %%%%%%%%%%%%%%%%%%%% x-position of colorbar   (+) --->>   +0.02;
        pos(1)=pos(1)+0.022; % 0.035
        pos(2)=pos(2)+0.09;  
        pos(3)=pos(3)+ 0.012;  % widthcolorbar
        pos(4)=pos(4)-0.18;  % height colorbar
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
        

        
        set(h,'Position',pos)

% if strcmp(iso,'dD')==1
%     iso_str='{\delta}D';
% elseif strcmp(iso,'Acc')==1
%     iso_str='Acc';
% elseif strcmp(iso,'MSA')==1
%     iso_str='MSA';
% elseif strcmp(iso,'d18O')==1    
%      iso_str='{\delta}^1^8O';
% elseif strcmp(iso,'d-excess')==1    
%      iso_str='d-excess';    
% end



   mlabel('MLabelParallel',-46  ,'PLabelLocation', [-75 -60 -45 -30],'PLabelMeridian',100) 


%%%%%%%%%%%%%%%%%%%%%%%%%%
% letter
TextBox = uicontrol('style','text');
% if  strcmp(proj,'merc')==1          
%             
%   letter_position=  [180 500 40 50];% [200 635 50 50];
%   
%     %set(TextBox,'String',letter,'position',[130 820 50 50],'FontWeight','bold','FontSize',28 ); % x position y position xsize ysize  
%     set(TextBox,'String',letter,'position',letter_position,'FontWeight','bold','FontSize',28 ); % x position y position xsize ysize            
%             
% elseif  strcmp(proj,'stereo')==1
    
       if strcmp(name,'SIC')==1 
       letter_position=  [200 820 60 60];%
       else
       %letter_position=   [330 840 50 50];
       letter_position=   [200 820 60 60];
       end
    
        
       
%        if yr_e==2009
            set(TextBox,'String',letter,'position',letter_position,'FontWeight','bold','FontSize',44 ); % x position y position xsize ysize 
%        else
%            set(TextBox,'String',letter,'position',letter_position,'FontWeight','bold','FontSize',40 ); % x position y position xsize ysize
%        end
            
            
            %285
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
  
    
  
  
%   if  season_nr==1
  % load(['C:\PHD\matlab_storage_of_output_files\Palmer_',iso,'_',z_level_str,'_',season,'_ERA_5_correlation _1979_2011stereo_crop.mat'])  
   load(['C:\PHD\matlab_storage_of_output_files\anomalies_2010_2011_z500_c2.mat'])   
%   else
      
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
    
     
     latg_c=latg((189:end),:); % only show winds inside figure limits
     long_c=long((189:end),:);
     Rc_cu_c=Rc_cu';
     Rc_cv_c=Rc_cv';
     
    wind_den=1;
     density_limit=10;
     unitString='ms^-^1';
    
     quivermc_c(double(latg_c), double(long_c), Rc_cu_c((189:end),:), Rc_cv_c((189:end),:),...
    'linewidth',1.5,'reference',wind_den,'density',density_limit,'units',unitString,'color',[.9 .0 .0]); % 'reference','median'
     
%  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 coast_nr=1;  % On/Off
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
%         if show_title==1
yr_s=2010;
yr_e=2011;            
%             if season_nr==1
            
                        hp1=  title(['Anomalies SIC, Z500, winds',...
                            ' ', num2str(yr_s),' - ',num2str(yr_e)],'FontSize',20, 'FontWeight', 'bold');
                        
%             elseif season_nr==6
%                                         hp1=  title([name,' ',SST_dataset_str,' ',season_str,...
%                             ' ', num2str(yr_s),' - ',num2str(yr_e)],'FontSize',18, 'FontWeight', 'bold');             
%                         
%             else
%                                         hp1=  title([name,' ',SST_dataset_str,' ',season_str,...
%                             ' ', num2str(yr_s+1),' - ',num2str(yr_e)],'FontSize',18, 'FontWeight', 'bold');
%             end
                        
                        
                        
%         elseif show_title==2
%                         hp1=  title([iso_str,', ',name],'FontSize',28, 'FontWeight', 'bold');                                                
%         elseif show_title==3
%                         hp1=  title([iso_str,' ',name,' ',SST_dataset_str,' ',...
%                             ' ', num2str(round(date_annual(start_t))),' - ',num2str(round(date_annual(end_t)))],'FontSize',20, 'FontWeight', 'bold');
%         elseif show_title==4
%                         hp1=  title(['r(',iso_str,', ',name,')'],'FontSize',28, 'FontWeight', 'bold');  
%         end

 % move title       
        if strcmp(proj,'stereo')==1  
            pos=get(hp1,'Position');
            pos(2)=pos(2)-0.03;
            set(hp1,'Position',pos,'FontSize',28, 'FontWeight', 'bold');    
%         elseif strcmp(proj,'mec')==1
%             pos=get(hp1,'Position');
%             pos(2)=pos(2)-0.02;
%         	set(hp1,'Position',pos,'FontSize',18, 'FontWeight', 'bold');    
        end
       
end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if show_colorbar==1
        set(h, 'FontSize',20, 'FontWeight', 'bold'); 
    end
    
    
    
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%  
% 
% if corr_w_std==0 && movstd_iso==0
%  
%  SIC_box=1; 
%  
% elseif corr_w_std==1 && movstd_iso==0 
%  SIC_box=2; 
% end
%  
%  if SIC_box==1    % eastern Weddell Sea
%         
% %         lon_b_c=[-26 -26 -40 -40 -26];
% %         lat_b_c=[ -75 -72 -72 -75 -75];
% %         lon_b_c=[-29 -29 -38 -38 -29];
% %         lat_b_c=[ -74.5 -72.5 -72.5 -74.5 -74.5];
% 
%         
%         lon_b_c=[-65 -65 -93 -93 -65]; % Bellinghausen lat line
%         lat_b_c=[ -65 -65 -65 -65 -65];
%         plotm(lat_b_c,lon_b_c,'-','LineWidth',4,'Color',[1 .2 .2])  % 1 .7 .1   
%         
%         
%  elseif SIC_box==2      
%              
% %         lon_b_c=[-31 -31 -37 -37 -31];
% %         lat_b_c=[ -74.3 -73.5 -73.5 -74.3 -74.3];
% %         
% % %         lon_b_c=[-31 -31 -36 -36 -31];
% % %         lat_b_c=[ -74.1 -73.5 -73.5 -74.1 -74.1];
% %         
% %         plotm(lat_b_c,lon_b_c,':','LineWidth',4,'Color',[1 .2 .2])       
% 
% 
%         lon_b_c=[-31 -31 -37 -37 -31];
%          lat_b_c=[ -74.3 -73.5 -73.5 -74.3 -74.3];
% %        lat_b_c=[ -74.2 -73.0 -73.0 -74.2 -74.2]; % not stronger corr
%         
%         plotm(lat_b_c,lon_b_c,':','LineWidth',4,'Color',[1 .7 .1])       
% end 



%%%%%%%%%%%%%%%%%%%%%%%%

%  if corr_w_std==1
%      
%     plotm(ndsic_lat(in_or1,in_or2),ndsic_lon(in_or1,in_or2),'o','color',color_id2,'MarkerSize',15,'LineWidth', 2);
%     
%     plotm(ndsic_lat(in_or3,in_or4),ndsic_lon(in_or3,in_or4),'o','color',color_id1,'MarkerSize',15,'LineWidth', 2);
%     
%  end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Combine Trends, SIC (shading), Z500 (contours), Winds (vectors)
% combined_fields_figure_Palmer_corr_SIC_z500_winds.m
% DE 2020
%
% MATLAB 2018a
% uses saved fields from:
% * Regression_spatial_ERA_I_z500_SAT_DE20_c17.m
% * Regression_spatial_HadISST_SIC_DE20_c6.m
% Github version 1   []
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
clear all
close all


corr_w_std=0; %(1/0) std variable, 0 mean 
movstd_iso=0; % (1/0) moving std isotopes
movstd_win=16; %11
%movstd_win=21;
type_nr=2;
period_nr=2;



% season_nr=2; % %(1) annual (2) DJF (3) MAM (4) JJA (5) SON

for season_nr=1:1 % change to 5 to do all seasons
    
    
    if period_nr==1
    yr_s=1979;
    yr_e=2011;
    step_nr=1;
    wind_den=0.2;
    elseif period_nr==2
    yr_s=2000;
    yr_e=2011;
    step_nr=5;
    wind_den=1;
    end

if season_nr==1
    season='annual';
    
elseif season_nr==2
     season='DJF';
     
elseif season_nr==3
     season='MAM';
     
elseif season_nr==4
     season='JJA';
     
elseif season_nr==5
     season='SON';
     
elseif season_nr==6
     season='MAMJJA';    
end


w_level_nr=3; % (1) 200 (2) 500 (3) 850       for V U winds
z_level_nr=1;

 iso_nr=4;
  
    if iso_nr==1
    iso='dD';
    elseif iso_nr==2
    iso='Acc';  
    elseif iso_nr==3
    iso='MSA';
    elseif iso_nr==4
    iso='d18O';
    elseif iso_nr==5
    iso='d-excess';    
    end

    name='SIC';
    if corr_w_std==1
    name=['std ',name];
    end   
    
    if movstd_iso==1
        iso_str_c=['mstd ',iso]; 
        iso_str_c2=['mstd_',num2str(movstd_win),'_',iso]; % for file name    
    else
        iso_str_c=iso;
        iso_str_c2=iso;
    end
        

   load(['C:\PHD\matlab_storage_of_output_files\HadISST_regression_SIC_Trend__monthly_2000-2011.mat']) 


tic
% Defaults
  proj='stereo';
  lat1=-90;
  
  
 %  lat2=-50;
  lat2=-50; 
%    lat2=-10; 
  
  
  box_use=0;
  
  
  if season_nr==1
      if period_nr==1
        letter='a';
      elseif period_nr==2
        letter='a';  
      end
  elseif season_nr==2
  letter='a';
  elseif season_nr==3
  letter='b';
  elseif season_nr==4
  letter='c';  
  elseif season_nr==5
  letter='d';
  elseif season_nr==6
      
   if corr_w_std==0   
    letter='a';
   elseif corr_w_std==1 && movstd_iso==1
    letter='c';  
   elseif corr_w_std==1
    letter='b';  
   end
  
  end
  
  
%   name='SIC';
 
  show_maximum_point=0;
  show_title=1;
  
  lock_scalebar=1;%%%%%%%%%%%%%%%%%%%%%%%
  
  Int_in=0;
  
%   if Int_in==1
%     RSAS_SIC_box=0;    
%   elseif Int_in==0    
%     RSAS_SIC_box=0; % (1) area defined using RICE corr pattern (2) ADP following Yuan 2004
%   end
  
  if season_nr==1
     season_str='Annual';
  elseif season_nr==2
     season_str='DJF';
  elseif season_nr==3
     season_str='MAM';
  elseif season_nr==4
     season_str='JJA';
  elseif season_nr==5
     season_str='SON';
  elseif season_nr==6
     season_str='MAMJJA';    
  end
  
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


Rc=s_all;

Rc5=s_all;
% Rc2=R;

% Rc(P>=p_level_rc)=NaN; 
% Rc2(P>=0.01)=NaN; % Mask out above or equal R's
Rc5(p_all>=0.1)=NaN; 


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
  X_r=squeeze(Rc(1,2,:,:))';
  X_r(:,(161:316))=NaN;
  [X_r_min, min_ind_1]=nanmin(nanmin(X_r)); % negative corr area RS
  [X_r_min_c, min_ind_2]=nanmin(nanmin(X_r,[],2));
  

  Px_c=squeeze(p_all(1,2,:,:))';
  
  p_level(Px_c(min_ind_2,min_ind_1)); % p level at correlation extreme

  
if lock_scalebar==1 
   
     c_limit=15; % 
         
end
 

Rc6=squeeze(Rc5(1,2,:,:))';
% Rc6(:,361)=Rc6(:,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if type_nr==1
    wr_c=s_all*10;
elseif  type_nr==2   
    wr_c=s_all*10*12*100; % per dec
end

hSurf=surfm(double(HadISST_lat),double(HadISST_lon),wr_c');
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

 toc

%%%%%%%%%%%%%%%%%%%%
% contourm function keeps projection defined above


contour_alt_nr=2; % (1), regular p (2) including Wilks test
if contour_alt_nr==1
    h1= contourm(double(HadISST_lat),double(HadISST_lon),p_all', [0.05],'-k','ShowText','off',...
   'LineWidth', 2);
elseif contour_alt_nr==2
    h1= contourm( double(HadISST_lat(146:167)),double(HadISST_lon),p_fdr_rej', [0.05],'-k','ShowText','off',...
    'LineWidth',2.5);
end


% Then reset the color axis according to the range determined above:
caxis(cRange); 
alpha(1) 
%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%


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


   % a1=axestext_c(1.09,+0.6, ['Regression'] ,'rotation',-90,'FontSize',22, 'FontWeight', 'bold');    %050
    a1=axestext_c(1.09,+0.6, ['Trend (% dec^-^1)'] ,'rotation',-90,'FontSize',22, 'FontWeight', 'bold'); 

      
    
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

if strcmp(iso,'dD')==1
    iso_str='{\delta}D';
elseif strcmp(iso,'Acc')==1
    iso_str='Acc';
elseif strcmp(iso,'MSA')==1
    iso_str='MSA';
elseif strcmp(iso,'d18O')==1    
     iso_str='{\delta}^1^8O';
elseif strcmp(iso,'d-excess')==1    
     iso_str='d-excess';    
end



   mlabel('MLabelParallel',-46  ,'PLabelLocation', [-75 -60 -45 -30],'PLabelMeridian',100) 


%%%%%%%%%%%%%%%%%%%%%%%%%%
% letter
TextBox = uicontrol('style','text');
    
       if strcmp(name,'SIC')==1 
       letter_position=  [200 820 60 60];%
       else
       %letter_position=   [330 840 50 50];
       letter_position=   [200 820 60 60];
       end
    
        
       
       if yr_e==2009
            set(TextBox,'String',letter,'position',letter_position,'FontWeight','bold','FontSize',44 ); % x position y position xsize ysize 
       else
           set(TextBox,'String',letter,'position',letter_position,'FontWeight','bold','FontSize',40 ); % x position y position xsize ysize
       end
            
            set(TextBox,'foregroundcolor', [0 0 0], ...
            'backgroundcolor', [1 1 1]);
               
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
  %clear all
    if z_level_nr==1
       z_level_str='z200';
    elseif z_level_nr==2
       z_level_str='z500';  
    elseif z_level_nr==3
       z_level_str='z850'; 
    end
  
      
   load(['C:\PHD\matlab_storage_of_output_files\ERA_interim_regression_z500_annual_Monthly_2000-2011.mat'])   


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
%   %contourm function keeps projection defined above
%     Pc=squeeze(P(1,2,:,:))';
if type_nr==1
    Rc_c1=s_all*10;
elseif type_nr==2
   Rc_c1=s_all*10*12; 
end

if type_nr==1
Rc_cc1=s_all*10;
elseif type_nr==2
Rc_cc1=s_all*10*12;    
end
Rc_c1(s_all>0)=NaN; % dashed

Rc_cc1(s_all>0)=NaN;% line contour
Rc_cc1(p_all>=0.05)=NaN;  % 


Pc_c_1=p_all;

Pc_c_1(s_all>0)=NaN;


%%%%%%
if type_nr==1
Rc_c2=s_all*10;
elseif type_nr==2
Rc_c2=s_all*10*12;    
end
Rc_c2(s_all<0)=NaN;


if type_nr==1
Rc_cc2=s_all*10;
elseif type_nr==2
Rc_cc2=s_all*10*12;    
end

Rc_cc2(s_all<0)=NaN;
Rc_cc2(p_all>=0.05)=NaN;


Pc_c_2=p_all;

Pc_c_2(s_all<0)=NaN;


     contourm( double(era_lat),double(era_long),Rc_c1','LevelStep',step_nr, 'LineColor',[0.0 0.2 0.6],'LineStyle',':','ShowText','off','LineWidth', 2); 
  
     if nanmin(nanmin(Pc_c_1)) <=0.05
     contourm( double(era_lat),double(era_long),Rc_cc1','LevelStep',step_nr, 'LineColor',[0.0 0.2 0.6],'LineStyle','-','ShowText','off','LineWidth', 2);
     end
     
     contourm( double(era_lat),double(era_long),Rc_c2','LevelStep',step_nr, 'LineColor',[1.0 0.1 0],'LineStyle',':','ShowText','off','LineWidth', 2);
 
% not totaly reliable (reason? sometimes to few point need more than one point to make contour), so comment/uncomment
% not working for annual
     if nanmin(nanmin(Pc_c_2)) <=0.05 %&& season_nr>1 % && ~(season_nr==6 && corr_w_std==0 && movstd_iso==0)
      contourm( double(era_lat),double(era_long),Rc_cc2','LevelStep',step_nr, 'LineColor',[1.0 0.1 0],'LineStyle','-','ShowText','off','LineWidth', 2);    
     end
     
     %%%%%%%%%%%%%%%%%%%%%
     if contour_alt_nr==2
                 co_pa=[.9 .1 .9];
                 line_w_nr=2;
         
            h1= contourm( double(era_lat),double(era_long),p_fdr_rej', [0.05],'-','ShowText','off',...
    'Linecolor',co_pa,'LineWidth', 4);

     end
     %%%%%%%%%%%
     
% hold on

show_wind_corr=1;

 % wind corr
 if show_wind_corr==1
       
    
    if w_level_nr==2
       w_level_str='v500';  
    elseif w_level_nr==3
       w_level_str='v850'; 
    end
    
    
    if corr_w_std==1
    w_level_str=['std ',w_level_str];
    end  
    
    
    load(['C:\PHD\matlab_storage_of_output_files\ERA_interim_regression_v850_annual_Monthly_2000-2011.mat'])
    
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    clear w_level_str
    
    if type_nr==1
    Rc_cv=s_all*10;
    elseif type_nr==2
    Rc_cv=s_all*10*12;    
    end 
    Pc_v=p_all;
   % Pc_v=p_fdr_rej;
    

    if w_level_nr==2
       w_level_str='u500';  
    elseif w_level_nr==3
       w_level_str='u850'; 
    end  
    
    if corr_w_std==1
    w_level_str=['std ',w_level_str];
    end
    
    
    clear  s_all
  
     load(['C:\PHD\matlab_storage_of_output_files\ERA_interim_regression_u850_annual_Monthly_2000-2011.mat'])
  
     
       %%%%%%%%%%%%%%%%%%%%%
     if contour_alt_nr==2
                 co_pa=[1.0 .7 .0];
                 line_w_nr=2;
         
            h1= contourm( double(era_lat),double(era_long),p_fdr_rej', [0.05],'-','ShowText','off',...
    'Linecolor',co_pa,'LineWidth', 4);

     end
     
     
  %%%%%%%%%%%
  
  if type_nr==1
    Rc_cu=s_all*10;
  elseif type_nr==2
    Rc_cu=s_all*10*12;  
  end
  
  
    Pc_u=p_all;
    
   % Pc_u=p_fdr_rej;
    
    %%%%%%%%%% Only show significant wind r
    
    [a b c]=size(Rc_cv);
    
    ind_w_v=zeros(a,b);
    ind_w_u=zeros(a,b);
%     ind_w_tot=NaN(a,b);
    
    for i=1:a
       for j=1:b
        
           if Pc_v(i,j)<0.05
               ind_w_v(i,j)=1;              
           end
           
           if Pc_u(i,j)<0.05
               ind_w_u(i,j)=1;              
           end          
                   
       end
    end
    
    
    ind_w_tot=ind_w_v+ind_w_u;
    
    Rc_cu(ind_w_tot==0)=NaN;
    Rc_cv(ind_w_tot==0)=NaN; 
    %     ind_v=find( Rc_cv(Pc>=0.1));
    %%%%%%%%
    
    clear  Rc_c
    
    
     [long,latg] = meshgrid(era_long,era_lat );
    
     
     latg_c=latg((52:90),:); % only show winds inside figure limits
     long_c=long((52:90),:);
     Rc_cu_c=Rc_cu';
     Rc_cv_c=Rc_cv';
    
     density_limit=15;
     unitString='ms^-^1/dec';
    
     quivermc_c(double(latg_c), double(long_c), Rc_cu_c((52:90),:), Rc_cv_c((52:90),:),...
    'linewidth',1.5,'reference',wind_den,'density',density_limit,'units',unitString,'color',[.9 .0 .0]); % 'reference','median'
     
 end
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% title
if show_title==1 || show_title==2 || show_title==3 || show_title==4
                                   % SST_dataset_str='ndsic';
                                   % SST_dataset_str='NSIDC';   
                                    SST_dataset_str='HadISST';   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
        if show_title==1
            
            if season_nr==1
            
                        hp1=  title(['Trend ' 'SIC, Z500, Winds',' ',season_str,...
                            ' ', num2str(yr_s),' - ',num2str(yr_e)],'FontSize',20, 'FontWeight', 'bold');
                        
            elseif season_nr==6
                                        hp1=  title([name,' ',SST_dataset_str,' ',season_str,...
                            ' ', num2str(yr_s),' - ',num2str(yr_e)],'FontSize',18, 'FontWeight', 'bold');             
                        
            else
                                        hp1=  title([name,' ',SST_dataset_str,' ',season_str,...
                            ' ', num2str(yr_s+1),' - ',num2str(yr_e)],'FontSize',18, 'FontWeight', 'bold');
            end
                        
                        
                        
        elseif show_title==2
                        hp1=  title([iso_str,', ',name],'FontSize',28, 'FontWeight', 'bold');                                                
        elseif show_title==3
                        hp1=  title([iso_str,' ',name,' ',SST_dataset_str,' ',...
                            ' ', num2str(round(date_annual(start_t))),' - ',num2str(round(date_annual(end_t)))],'FontSize',20, 'FontWeight', 'bold');
        elseif show_title==4
                        hp1=  title(['r(',iso_str,', ',name,')'],'FontSize',28, 'FontWeight', 'bold');  
        end

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

 if corr_w_std==1
     
    plotm(ndsic_lat(in_or1,in_or2),ndsic_lon(in_or1,in_or2),'o','color',color_id2,'MarkerSize',15,'LineWidth', 2);
    
    plotm(ndsic_lat(in_or3,in_or4),ndsic_lon(in_or3,in_or4),'o','color',color_id1,'MarkerSize',15,'LineWidth', 2);
    
 end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

hold off
%%%%%%%%%%%%%%%%%%%%
% save figure
 filedir ='C:\PHD\matlab_storage_of_output_files\figures\'; 
                filename=strcat(filedir,'combined_trend_HadISST_ERAI_',name,'_',season_str,'_',z_level_str,'_',...
                    w_level_str,'_',proj,'_',num2str(yr_s),'_',num2str(yr_e));
                
savefilename_c=filename;

 % save as png 
orient landscape
% increase quality of picture by increase to -r500
% if figure_format==1
%   export_fig('-eps','-nocrop','-painters', '-depsc','-opengl', '-r100',savefilename_c); % EPS works  func needs ghostscript 
% elseif figure_format==2
%  export_fig('-png','-painters', '-depsc','-nocrop','-opengl', '-r190',savefilename_c); % PNG  110
  export_fig('-pdf','-painters', '-depsc','-nocrop','-opengl', '-r190',savefilename_c); % PNG  110 
% end

end

toc
%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%
% SAT_SIC_z500_RICE_DE20_c10.m
% d18O corr w. SAT, SIC and z500
% Daniel Emanuelsson, 2020
% MATLAB 2018a 
% Github version 2 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions
% * IPO_vline.m DE
% * fig.m fileexchange
% * export_fig.m fileexchange
%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
%%%%%
addpath C:\PHD\matlab_lib\Data\
addpath C:\PHD\ERA_interim\
%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  
%     %load('C:\PHD\matlab_storage_of_output_files\RICE_combined_Deep_1213B_c19.mat'); % May 2017
    load('C:\PHD\matlab_storage_of_output_files\RICE_combined_Deep_1213B_c23.mat');
    stacked_record_annual_Ma=MA_save((1901:end),:);
    %iso_str='stack';
    iso_str='d18O';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fig   %%% moving correlation with RICE dD
sliding_w=11;%%%%%%%%%%%%%%%%%%%%%%%%%%%
add_yr=2; % (0) 2009 (2) 2011
    

h=fig('units','inches','width',14,'height',4.5,'font','Helvetica','fontsize',16,'border','on');
   box on
   hold on   
   set(gca,'YLim',[-1 1])
 %%%%%%%%%%%%%%            
addpath C:\PHD\matlab_lib\Library
%   edit IPO_vline
  IPO_vline(15)
   hold on 
   
load('C:\PHD\matlab_storage_of_output_files\Ma_ADP_index_HadISST_SIC.mat'); % ADP index   ADP_index.m  
  
load('C:\PHD\matlab_storage_of_output_files\SST_Nino-4_HadISST_1979_2012.mat');
%   ENSO 
  
  
      %%%% ERA-Interim
      %%%%%%%%%%%%%%%%%%%%%%% RICE
le=length(stacked_record_annual_Ma(80:110+add_yr));

Ma_corr2=NaN(le,3); % use NaN instead of zero because NPGO index is shorter

% iso
Ma_corr2(:,1)=detrend(stacked_record_annual_Ma((80:110+add_yr),3));
Ma_corr2(:,2)=detrend(Ma_save_nino(1:31+add_yr));

Ma_corr2(:,3)=detrend(MA_ADP_save((8:38+add_yr),2));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% variable 1: Nino-4 SST
% not using function to get better overview
% 

windowSize=sliding_w;

[N,M] = size(Ma_corr2);

correlationTS = nan(N, 1);
correlationTS_p= nan(N, 1);
correlationTS_rlo=nan(N, 1);
correlationTS_rup=nan(N, 1);


% correlationTS = nan(N, M-1);
% correlationTS_p= nan(N, M-1);
% correlationTS_rlo=nan(N, M-1);
% correlationTS_rup=nan(N, M-1);


indexColumn=1;
 %for t = windowSize+1:N % (Original)
for t = floor(windowSize/2)+1:(N-floor(windowSize/2)) % Modified
        % 6:105:
     
%      dataMatrix=Ma_corr4;
 %C = corrcoef(dataMatrix(t-windowSize:t, :)); %  (Original)
 % C = corrcoef_df(dataMatrix(t-windowSize:t, :)); % switch to UW Steig corr function
 % [C, p, rlo, rup ]= corrcoef_df(dataMatrix(t-windowSize:t, :)); % switch to UW Steig corr function 
 
 % [C, p, rlo, rup ]= corrcoef_df(dataMatrix(t-floor(windowSize/2):t+floor(windowSize/2), :)); % Change to centered moving correlation window
  %[C, p, rlo, rup ]= corrcoef(dataMatrix(t-floor(windowSize/2):t+floor(windowSize/2), :), 'alpha', alpha); 
  
  X=Ma_corr2(t-floor(windowSize/2):t+floor(windowSize/2),1);
  Y=Ma_corr2(t-floor(windowSize/2):t+floor(windowSize/2),2); 
  
  
  [C, p, rlo, rup ]= corrcoef_df(X,Y); % Change to centered moving correlation window 
 
 
%     idx = setdiff(1:M, [indexColumn]);
%     correlationTS(t, :) = C(indexColumn, idx);
%     correlationTS_p(t, :)= p(indexColumn, idx);
%     correlationTS_rlo(t, :)= rlo(indexColumn, idx);
%     correlationTS_rup(t, :)= rup(indexColumn, idx);
    
    
%    idx = setdiff(1:M, [indexColumn]);
    correlationTS(t, :) = C(1, 2);
    correlationTS_p(t, :)= p(1, 2);
    correlationTS_rlo(t, :)= rlo(1, 2);
    correlationTS_rup(t, :)= rup(1, 2);
     
    
end


clear X Y t C p rlo rup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% variable 2: ADP



correlationTS_c = nan(N, 1);
correlationTS_p_c= nan(N, 1);
correlationTS_rlo_c=nan(N, 1);
correlationTS_rup_c=nan(N, 1);


for t = floor(windowSize/2)+1:(N-floor(windowSize/2)) % Modified
        % 6:105:
     
%      dataMatrix=Ma_corr4;
 %C = corrcoef(dataMatrix(t-windowSize:t, :)); %  (Original)
 % C = corrcoef_df(dataMatrix(t-windowSize:t, :)); % switch to UW Steig corr function
 % [C, p, rlo, rup ]= corrcoef_df(dataMatrix(t-windowSize:t, :)); % switch to UW Steig corr function 
 
 % [C, p, rlo, rup ]= corrcoef_df(dataMatrix(t-floor(windowSize/2):t+floor(windowSize/2), :)); % Change to centered moving correlation window
  %[C, p, rlo, rup ]= corrcoef(dataMatrix(t-floor(windowSize/2):t+floor(windowSize/2), :), 'alpha', alpha); 
  
  X=Ma_corr2(t-floor(windowSize/2):t+floor(windowSize/2),1);
  Y=Ma_corr2(t-floor(windowSize/2):t+floor(windowSize/2),3); 
  
  
  [C, p, rlo, rup ]= corrcoef_df(X,Y); % Change to centered moving correlation window 
 
 
%     idx = setdiff(1:M, [indexColumn]);
%     correlationTS(t, :) = C(indexColumn, idx);
%     correlationTS_p(t, :)= p(indexColumn, idx);
%     correlationTS_rlo(t, :)= rlo(indexColumn, idx);
%     correlationTS_rup(t, :)= rup(indexColumn, idx);
    
    
%    idx = setdiff(1:M, [indexColumn]);
    correlationTS_c(t, :) = C(1, 2);
    correlationTS_p_c(t, :)= p(1, 2);
    correlationTS_rlo_c(t, :)= rlo(1, 2);
    correlationTS_rup_c(t, :)= rup(1, 2);
     
    
end


clear X Y t C p rlo rup

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  [correlationTS2, correlationTS_p2, correlationTS_rlo2, correlationTS_rup2]=movingCorrelation_c(Ma_corr2,sliding_w,1, 0.05);    


% edit corrcoef
% edit corrcoef_df
% edit movingCorrelation_c


color_code=[0.3,0.3,0.3];

h70=plot(stacked_record_annual_Ma((80:110+add_yr),1), correlationTS_rup(:,1),'--','LineWidth',2,'MarkerSize',14,'Color',color_code);
h7=plot(stacked_record_annual_Ma((80:110+add_yr),1), correlationTS(:,1),'-','LineWidth',4,'MarkerSize',14,'Color',color_code);
h71=plot(stacked_record_annual_Ma((80:110+add_yr),1), correlationTS_rlo(:,1),'--','LineWidth',2,'MarkerSize',14,'Color',color_code);

% h80=plot(stacked_record_annual_Ma((80:110+add_yr),1), correlationTS_rup2(:,2)*(-1),'--','LineWidth',2,'MarkerSize',14,'Color',[0.6,0.0,0.0]); 
% h8=plot(stacked_record_annual_Ma((80:110+add_yr),1), correlationTS2(:,2)*(-1),'-','LineWidth',4,'MarkerSize',14,'Color',[0.6,0.0,0.0]);
% h81=plot(stacked_record_annual_Ma((80:110+add_yr),1), correlationTS_rlo2(:,2)*(-1),'--','LineWidth',2,'MarkerSize',14,'Color',[0.6,0.0,0.0]); 


h80=plot(stacked_record_annual_Ma((80:110+add_yr),1), correlationTS_rup_c(:,1)*(-1),'--','LineWidth',2,'MarkerSize',14,'Color',[0.6,0.0,0.0]); 
h8=plot(stacked_record_annual_Ma((80:110+add_yr),1), correlationTS_c(:,1)*(-1),'-','LineWidth',4,'MarkerSize',14,'Color',[0.6,0.0,0.0]);
h81=plot(stacked_record_annual_Ma((80:110+add_yr),1), correlationTS_rlo_c(:,1)*(-1),'--','LineWidth',2,'MarkerSize',14,'Color',[0.6,0.0,0.0]); 

set(gca,...
          'linewidth',3,...
          'FontWeight','bold' );
    
   %%%%% Legend %%%%%%%%%%%%
    hs61='PC1 SAT';
    hs89='PC2 SAT';
    hs62='PC3 SAT';
    hs7='Nino-4 SST'; 
    hs8='ADP (-1)'; 
        
    hs7_c=['r ({\delta}^1^8O, ',hs7,')'];
    hs8_c=['r ({\delta}^1^8O, ',hs8,')'];

  h_leg=legend([ h7 h8 ],  hs7_c, hs8_c);
 %  h_leg=legend([  h8 ],   hs8_c);
  set(h_leg, 'location', 'SouthWest','EdgeColor',[0 0 0],'FontSize',18,'FontWeight','bold' ); % ,'color','none'
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % letter
  letter_pos=[80 350 60 60];
  letter='d';
          TextBox = uicontrol('style','text');
          set(TextBox,'String',letter,'position',letter_pos,'FontWeight','bold','FontSize',34 ); % x position y position xsize ysize
          set(TextBox,'foregroundcolor', [0 0 0], ...
         'backgroundcolor', [1 1 1]);
 %%%%%%%%%%%%%%%%%%%%%%%%%
        set(gca,'XLim',[1977 2012])
        %hl=hline([-0.5,0,0.5]);
        hl=hline([0]);
        set(hl,'LineWidth',3);
        
%         xlabel('Year','FontWeight','bold','FontSize',18 );
         ylabel('Correlation','FontWeight','bold','FontSize',20 );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
% save fig.
end_str='_c2';
filename=['ADP_RICE_w_',num2str(sliding_w),'_', iso_str, end_str]; % poor chooise of name
filedir ='C:\PHD\matlab_storage_of_output_files\figures\';
savefilename_c=strcat(filedir,filename);
% save as png 
orient landscape
% increase quality of picture by increase to -r300
%export_fig('-png','-nocrop','-painters', '-depsc','-opengl', '-r170',savefilename_c); % PNG-nocrop' 
export_fig('-pdf','-nocrop','-painters', '-depsc','-opengl', '-r190',savefilename_c); % PNG-nocrop'

%% scatter plots
close all


%alt=5; % (1) ADP (2) nino-4 (3) SAM (4) PSA1 (5) PSA2

for alt=1:5

 h=fig('units','inches','width',6,'height',6,'font','Helvetica','fontsize',16,'border','on');  
   box on
   hold on
   
  x1=stacked_record_annual_Ma((80:110+add_yr),3);
  
  if alt==1
    y1=MA_ADP_save((8:38+add_yr),2);
  elseif alt==2
    y1=Ma_save_nino(1:31+add_yr);
  elseif alt>2
      
      load('C:\PHD\matlab_storage_of_output_files\ERA-Interim_PCs_z500_lim0-360_-20_-90_1979-2012_annual_c5_annual_mean_rotate.mat') %z500
    
      if alt==3
        pc_nr=2;
        factr=-1;
      elseif alt==4
        pc_nr=3;
        factr=-1;
      elseif alt==5
        pc_nr=4;
        factr=-1;
      end
      y1=MA_PCs_save((1:33), pc_nr)*factr;
      
      if alt==2
          
      end
    
  end
   
plot(x1,y1,'*','MarkerSize',10,'Color',[0.0,0.0,0.0])

 if alt==1 || alt==2
    plot(x1(end:end),y1(end:end),'*','MarkerSize',10,'Color',[1,0.0,0.0])
 elseif alt==4
    plot(x1(end-1:end),y1(end-1:end),'*','MarkerSize',10,'Color',[1,0.0,0.0]) 
 end


           set(gca,...
          'linewidth',3,...
          'FontWeight','bold' );
xlabel(['{\delta}^1^8O (',char(8240),')'])



if alt==1 || alt==2  
x1=x1(1:end-1);
y1=y1(1:end-1);

elseif alt==4
    
x1=x1(1:end-2);
y1=y1(1:end-2);    

end

pol_1=polyfit(x1,y1,1);
x2=[-29:0.5:-15];
polv_1=polyval(pol_1,x2);

plot(x2,polv_1,'-k','LineWidth',2)
set(gca,'XLim',[-28 -18])

if alt==1
    
    
  title_str='ADP SIC';
    
elseif alt==2
    
set(gca,'YLim',[-0.8 0.8])  

 title_str='Nino-4 SST (°C)';
    
elseif alt==3 
  title_str='SAM (m s.d.^-^1)';  
elseif alt==4 
  title_str='PSA1 (m s.d.^-^1)';
elseif alt==5 
  title_str='PSA2 (m s.d.^-^1)';  
 
end

ylabel(title_str)


%correlation
%  factor_nr
[r,p, nv]=corrcoef_df_c2(detrend(x1),detrend(y1));
rp=[r(2),p(2)]
p_level(p(2))
%nv
% edit corrcoef_df_c2
% help getdof

%  x=X;
%  y=Y;

 nc=size(x1);
   r1 = corrcoef(x1(2:end),x1(1:end-1),'rows','pairwise');
   r2 = corrcoef(y1(2:end),y1(1:end-1),'rows','pairwise');
%    r1 = corrcoef(x(2:end),x(1:end-1),'rows','pairwise');
%    r2 = corrcoef(y(2:end),y(1:end-1),'rows','pairwise');
%    %r1 = acf(x,1);
%    %r2 = acf(y,2);
%    
  %  nv = nv*(1-r1(2)*r2(2))/(1+r1(2)*r2(2))
  
    neff = round(nc(1)*(1-r1(2)*r2(2))/(1+r1(2)*r2(2)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if alt==1
x_text=0.1;    
y_text=0.1;
letter='e';
elseif alt==2
y_text=0.98;
x_text=0.1;
letter='f';
elseif alt==3
y_text=0.98;
x_text=0.9;
letter='b';
elseif alt==4
y_text=0.98;
x_text=0.1;
letter='c';
elseif alt==5
y_text=0.98;
x_text=0.1;
letter='d';
end



 axestext_c(x_text,y_text,['r = ',num2str(round(rp(1)*100)/100),', ',p_level(p(2)),', n_e_f_f = ',num2str(neff)],'FontWeight','bold','FontSize',14 );   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % letter
  letter_pos=[0 485 35 65];
%   letter='b';
          TextBox = uicontrol('style','text');
          set(TextBox,'String',letter,'position',letter_pos,'FontWeight','bold','FontSize',30 ); % x position y position xsize ysize
          set(TextBox,'foregroundcolor', [0 0 0], ...
         'backgroundcolor', [1 1 1]);
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%          
% save fig.         
filename=['scatter_plot_',num2str(alt),'_', iso_str,'_scatter_plot']; % poor chooise of name
filedir ='C:\PHD\matlab_storage_of_output_files\figures\';
savefilename_c=strcat(filedir,filename);
% save as png 
orient landscape
% increase quality of picture by increase to -r300
%export_fig('-png','-nocrop','-painters', '-depsc','-opengl', '-r170',savefilename_c); % PNG-nocrop' 
export_fig('-pdf','-painters', '-depsc','-opengl', '-r190',savefilename_c); % PNG-nocrop'

end

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Z500
% 
    color_code_1=[0.1,0.1,0.1];
    color_code_2= [0.9,0.0,0.0];
    color_code_3= [0.0,0.0,0.8];

load('C:\PHD\matlab_storage_of_output_files\ERA-Interim_PCs_z500_lim0-360_-20_-90_1979-2012_annual_c5_annual_mean_rotate.mat') %z500

MA_PCs_save_z500=MA_PCs_save;  
% MA_PCs_save_monthly_z500=MA_PCs_save;
clear MA_PCs_save MA_PCs_save_monthly  
%%%%%%%%%%%%%%%%%%%%%%% RICE
% pol_nr11=1; % polarity positive SAM, as displayed in Fig 2 2014
% pol_nr12=1; % PSA patterns as in Kidson 1988 his fig 4b,c 
% pol_nr13=-1;
pol_nr11=-1; % polarity positive SAM, as displayed in Fig 2  2012
pol_nr12=-1; % PSA patterns as in Kidson 1988 his fig 4b,c 
pol_nr13=-1;

fac_nr11=1; 
fac_nr12=1;
fac_nr13=1;

le=length(stacked_record_annual_Ma(80:110+add_yr)); % 1979-2011
Ma_corr=NaN(le,6);
Ma_corr4(:,1)=detrend(stacked_record_annual_Ma((80:110+add_yr),2));
Ma_corr4((1:31+add_yr),2)=detrend(MA_PCs_save_z500((1:31+add_yr),2))*pol_nr11*fac_nr11; % z500 PC1 SAM related  (factor from Regression era)
Ma_corr4((1:31+add_yr),3)=detrend(MA_PCs_save_z500((1:31+add_yr),3))*pol_nr12*fac_nr12; % z500 PC 3 PSA2 related  
Ma_corr4((1:31+add_yr),4)=detrend(MA_PCs_save_z500((1:31+add_yr),4))*pol_nr13*fac_nr13; % z500 PC 3 PSA1 related

%edit movingCorrelation_c

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% variable 1: Z500 PC1
% not using function to get better overview
% 


clear correlationTS correlationTS_p correlationTS_rlo correlationTS_rup
clear correlationTS_c correlationTS_p_c correlationTS_rlo_c correlationTS_rup_c

windowSize=sliding_w;

[N,M] = size(Ma_corr4);

correlationTS = nan(N, 1);
correlationTS_p= nan(N, 1);
correlationTS_rlo=nan(N, 1);
correlationTS_rup=nan(N, 1);


% correlationTS = nan(N, M-1);
% correlationTS_p= nan(N, M-1);
% correlationTS_rlo=nan(N, M-1);
% correlationTS_rup=nan(N, M-1);


indexColumn=1;
 %for t = windowSize+1:N % (Original)
for t = floor(windowSize/2)+1:(N-floor(windowSize/2)) % Modified
        % 6:105:
     
%      dataMatrix=Ma_corr4;
 %C = corrcoef(dataMatrix(t-windowSize:t, :)); %  (Original)
 % C = corrcoef_df(dataMatrix(t-windowSize:t, :)); % switch to UW Steig corr function
 % [C, p, rlo, rup ]= corrcoef_df(dataMatrix(t-windowSize:t, :)); % switch to UW Steig corr function 
 
 % [C, p, rlo, rup ]= corrcoef_df(dataMatrix(t-floor(windowSize/2):t+floor(windowSize/2), :)); % Change to centered moving correlation window
  %[C, p, rlo, rup ]= corrcoef(dataMatrix(t-floor(windowSize/2):t+floor(windowSize/2), :), 'alpha', alpha); 
  
  X=Ma_corr4(t-floor(windowSize/2):t+floor(windowSize/2),1);
  Y=Ma_corr4(t-floor(windowSize/2):t+floor(windowSize/2),2); 
  
  
  [C, p, rlo, rup ]= corrcoef_df(X,Y); % Change to centered moving correlation window 
 
 
%     idx = setdiff(1:M, [indexColumn]);
%     correlationTS(t, :) = C(indexColumn, idx);
%     correlationTS_p(t, :)= p(indexColumn, idx);
%     correlationTS_rlo(t, :)= rlo(indexColumn, idx);
%     correlationTS_rup(t, :)= rup(indexColumn, idx);
    
    
%    idx = setdiff(1:M, [indexColumn]);
    correlationTS(t, :) = C(1, 2);
    correlationTS_p(t, :)= p(1, 2);
    correlationTS_rlo(t, :)= rlo(1, 2);
    correlationTS_rup(t, :)= rup(1, 2);
     
    
end


clear X Y t C p rlo rup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% variable 2: Z500 PC2



correlationTS_c = nan(N, 1);
correlationTS_p_c= nan(N, 1);
correlationTS_rlo_c=nan(N, 1);
correlationTS_rup_c=nan(N, 1);


for t = floor(windowSize/2)+1:(N-floor(windowSize/2)) % Modified
        % 6:105:
     
%      dataMatrix=Ma_corr4;
 %C = corrcoef(dataMatrix(t-windowSize:t, :)); %  (Original)
 % C = corrcoef_df(dataMatrix(t-windowSize:t, :)); % switch to UW Steig corr function
 % [C, p, rlo, rup ]= corrcoef_df(dataMatrix(t-windowSize:t, :)); % switch to UW Steig corr function 
 
 % [C, p, rlo, rup ]= corrcoef_df(dataMatrix(t-floor(windowSize/2):t+floor(windowSize/2), :)); % Change to centered moving correlation window
  %[C, p, rlo, rup ]= corrcoef(dataMatrix(t-floor(windowSize/2):t+floor(windowSize/2), :), 'alpha', alpha); 
  
  X=Ma_corr4(t-floor(windowSize/2):t+floor(windowSize/2),1);
  Y=Ma_corr4(t-floor(windowSize/2):t+floor(windowSize/2),3); 
  
  
  [C, p, rlo, rup ]= corrcoef_df(X,Y); % Change to centered moving correlation window 
 
    
%    idx = setdiff(1:M, [indexColumn]);
    correlationTS_c(t, :) = C(1, 2);
    correlationTS_p_c(t, :)= p(1, 2);
    correlationTS_rlo_c(t, :)= rlo(1, 2);
    correlationTS_rup_c(t, :)= rup(1, 2);
     
    
end



clear X Y t C p rlo rup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% variable 3: Z500 PC3


correlationTS_c3 = nan(N, 1);
correlationTS_p_c3= nan(N, 1);
correlationTS_rlo_c3=nan(N, 1);
correlationTS_rup_c3=nan(N, 1);


for t = floor(windowSize/2)+1:(N-floor(windowSize/2)) % Modified
        % 6:105:
     
%      dataMatrix=Ma_corr4;
 %C = corrcoef(dataMatrix(t-windowSize:t, :)); %  (Original)
 % C = corrcoef_df(dataMatrix(t-windowSize:t, :)); % switch to UW Steig corr function
 % [C, p, rlo, rup ]= corrcoef_df(dataMatrix(t-windowSize:t, :)); % switch to UW Steig corr function 
 
 % [C, p, rlo, rup ]= corrcoef_df(dataMatrix(t-floor(windowSize/2):t+floor(windowSize/2), :)); % Change to centered moving correlation window
  %[C, p, rlo, rup ]= corrcoef(dataMatrix(t-floor(windowSize/2):t+floor(windowSize/2), :), 'alpha', alpha); 
  
  X=Ma_corr4(t-floor(windowSize/2):t+floor(windowSize/2),1);
  Y=Ma_corr4(t-floor(windowSize/2):t+floor(windowSize/2),4); 
  
  
  [C, p, rlo, rup ]= corrcoef_df(X,Y); % Change to centered moving correlation window 
 
    
%    idx = setdiff(1:M, [indexColumn]);
    correlationTS_c3(t, :) = C(1, 2);
    correlationTS_p_c3(t, :)= p(1, 2);
    correlationTS_rlo_c3(t, :)= rlo(1, 2);
    correlationTS_rup_c3(t, :)= rup(1, 2);
     
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  [correlationTS4, correlationTS4_p, correlationTS4_rlo, correlationTS4_rup]=movingCorrelation_c(Ma_corr4,sliding_w,1, 0.05);
 
 h=fig('units','inches','width',14,'height',4.5,'font','Helvetica','fontsize',16,'border','on');
   box on
   hold on 
   set(gca,'YLim',[-1 1])
 %%%%%%%%%%%%%%
shade_nr=1;
        if shade_nr==1
            
            addpath C:\PHD\matlab_lib\Library
            IPO_vline(10)
            hold on
        
        end
     % PC1
       plot(stacked_record_annual_Ma((80:110+add_yr),1), correlationTS_rlo(:,1),'--','LineWidth',2,'MarkerSize',14,'Color',color_code_1);
    h1=plot(stacked_record_annual_Ma((80:110+add_yr),1), correlationTS(:,1),'-','LineWidth',4,'MarkerSize',14,'Color',color_code_1);
       plot(stacked_record_annual_Ma((80:110+add_yr),1), correlationTS_rup(:,1),'--','LineWidth',2,'MarkerSize',14,'Color',color_code_1);
    
    % PC2
       plot(stacked_record_annual_Ma((80:110+add_yr),1), correlationTS_rlo_c(:,1),'--','LineWidth',2,'MarkerSize',14,'Color',color_code_2);
    h2=plot(stacked_record_annual_Ma((80:110+add_yr),1), correlationTS_c(:,1),'-','LineWidth',4,'MarkerSize',14,'Color',color_code_2);     
       plot(stacked_record_annual_Ma((80:110+add_yr),1), correlationTS_rup_c(:,1),'--','LineWidth',2,'MarkerSize',14,'Color',color_code_2);
       
    % PC3   
%        plot(stacked_record_annual_Ma((80:110+add_yr),1), correlationTS_rlo_c3(:,1),'--','LineWidth',2,'MarkerSize',14,'Color',color_code_3);
%     h3=plot(stacked_record_annual_Ma((80:110+add_yr),1), correlationTS_c3(:,1),'-','LineWidth',4,'MarkerSize',14,'Color',color_code_3);
%        plot(stacked_record_annual_Ma((80:110+add_yr),1), correlationTS_rup_c3(:,1),'--','LineWidth',2,'MarkerSize',14,'Color',color_code_3);
    
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
%    hs3='PSA2';
    
%     hs1_c=['r (dD, ',hs1,' (',num2str(fac_nr11),'))'];
%     hs2_c=['r (dD, ',hs2,')'];
%     hs3_c=['r (dD, ',hs3,')'];
    
    hs1_c=['r ({\delta}^1^8O, ',hs1,')'];% (',num2str(fac_nr11),'))'];
    hs2_c=['r ({\delta}^1^8O, ',hs2,')'];
 %   hs3_c=['r ({\delta}^1^8O, ',hs3,')'];
 
    
  %h_leg=legend([ h1 h2 h3],  hs1_c, hs2_c, hs3_c); 
  h_leg=legend([ h1 h2],  hs1_c, hs2_c); 
  set(h_leg, 'location', 'SouthWest','EdgeColor',[0 0 0],'FontSize',18,'FontWeight','bold' ); % ,'color','none'
  
   title_str='Z500 PCs';   
 axestext_c(0.9999999,1.0,title_str,'FontWeight','bold','FontSize',20 ); 
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % letter
  letter_pos=[80 350 60 60];
  letter='a';
          TextBox = uicontrol('style','text');
          set(TextBox,'String',letter,'position',letter_pos,'FontWeight','bold','FontSize',34 ); % x position y position xsize ysize
          set(TextBox,'foregroundcolor', [0 0 0], ...
         'backgroundcolor', [1 1 1]);
 %%%%%%%%%%%%%%%%%%%%%%%%%
        set(gca,'XLim',[1977 2012])
      %  hl=hline([-0.5,0,0.5]);
        hl=hline([0]);
        set(hl,'LineWidth',3);
%        xlabel('Year','FontWeight','bold','FontSize',18 );
        ylabel('Correlation','FontWeight','bold','FontSize',20 );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save fig
end_str='_c2';
filename=['z500_RICE_ERA_interim_mov_corr_w_',num2str(sliding_w),'_', iso_str,end_str];       
filedir ='C:\PHD\matlab_storage_of_output_files\figures\';
savefilename_c=strcat(filedir,filename);
% save as png 
orient landscape
% increase quality of picture by increase to -r300
% export_fig('-png','-nocrop','-painters', '-depsc','-opengl', '-r170',savefilename_c); % PNG-nocrop' 
export_fig('-pdf','-nocrop','-painters', '-depsc','-opengl', '-r190',savefilename_c);
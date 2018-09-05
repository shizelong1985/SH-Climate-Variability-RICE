%%%%%%%%%%%%
% stepwise_MLR_Margaret_RICE_AWS_monthly_DE18_c2.m
% Matlab 2017a
% Daniel Emanuelsson 2018
% Github version 1  04-09-2018
%
% subfunctions

% * annave.m          UoW online archive Atmospheric Science
% * corrcoef_df.m     UoW Steig
% * fig.m             fileexchange
% * export_fig.m      fileexchange
%%%%%%%%%%%%
%
% file for processing of AWS data: 
% * AWS_Margaret_10min_data_c2.m
% * AWS_Roosevelt_island_c5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all
close all

%% month-to-month corr


load('C:\PHD\matlab_storage_of_output_files\Margaret_monthly_SAT_anomaly_Ma.mat') % AWS_Margaret climoloty removed
% Margaret (Lazzara et al. 2012) 
% University of Madison
load('C:\PHD\matlab_storage_of_output_files\ERA-Interim_PCs_2mT_lim0-360_-30_-90_1979-2014_annual.mat');
name_str='SAT';
       
       
  alt_runs=1;%%%%%%%%%%%%%%%%%%%
  
  if alt_runs==1
       
    a_c=361;
    b_c=432;
    
    d_c=1; % Margaret
    e_c=72;
    
  elseif alt_runs==2 % Same range as RICE AWS record
  
    a_c=385;
    b_c=418;
      
    d_c=25; % Margaret
    e_c=58;
    
  end

%
y=Margaret_monthly_SAT_anomaly_Ma((d_c:e_c),6); % (:,6) clim already removed

MA_PCs_save_monthly_c=MA_PCs_save_monthly((a_c:b_c),2);


[X_PC,clim]=annave(MA_PCs_save_monthly_c);

[r,p]=corrcoef_df(detrend(y),X_PC);

rp=[r(2),p(2)]

p_level(p(2))
 
      
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% stepwise MLR and SAT
  
  [X_PC1,clim1]=annave(MA_PCs_save_monthly((a_c:b_c),2)); % remove clim from PCs
  [X_PC2,clim2]=annave(MA_PCs_save_monthly((a_c:b_c),3));
  [X_PC3,clim3]=annave(MA_PCs_save_monthly((a_c:b_c),4));
  
  X_PC_step=zeros(length(X_PC1(1:end)),3);
  X_PC_step(:,1)=X_PC1(1:end);
  X_PC_step(:,2)=X_PC2(1:end);
  X_PC_step(:,3)=X_PC3(1:end);
 
LM=stepwiselm(X_PC_step,y(1:end),'Upper','linear'); % PCs
  
cur_fitted=LM.Fitted; % Fitted (predicted) response values based on input data

% plotDiagnostics(LM)

Coeff_tabl=LM.Coefficients;
Coeff_tabl_c=table2array(Coeff_tabl);
Coeff_Names=LM.CoefficientNames;      
      
 X_PC1_cp=X_PC_step(:,1);

b=regress(y,X_PC1_cp);
ry2=corrcoef_df(y,X_PC1_cp*b);

% ry2(2)^2

name_str='AWS SAT';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Margaret AWS
% Figure
% Fig. 14a

h=fig('units','inches','width',14,'height',5,'font','Helvetica','fontsize',16,'border','on');

    h1=plot(MA_PCs_save_monthly((a_c:b_c),1),y,'-','LineWidth',3,'MarkerSize',8,'Color',[0.1,0.1,0.1]); 
      hold on 

     
      
   %%%%%%%%%%%%%%%%%
   set(gca,'XLim',[2008.9 2015])
   set(gca,'YLim',[-15 15])      
      
    line1=cur_fitted;

    h2=plot(MA_PCs_save_monthly((a_c:b_c),1),line1,'--*','LineWidth',3,'MarkerSize',8,'Color',[0.65,0.05,0.0]);
    
    h3=plot(MA_PCs_save_monthly((a_c:b_c),1),X_PC1_cp*b,'--','LineWidth',3,'MarkerSize',8,'Color',[0.95,0.45,0.1]);
    
   
    ry1=corrcoef_df(y,line1);
    
    % ry1(2)^2
    
 %%%%%%%%%%%%%
 % legend
 hs1_c='SAT Margaret';
 % PCs
 Coeff_tabl_c_str1= num2str(round(Coeff_tabl_c(1,1)*1000)/1000);
 Coeff_tabl_c_str2= num2str(round(Coeff_tabl_c(2,1)*1000)/1000);
 Coeff_tabl_c_str3= num2str(round(Coeff_tabl_c(3,1)*1000)/1000); 
 
%  hs2_c= strcat(" ",Coeff_tabl_c_str2," ", Coeff_Names(2),' ',Coeff_tabl_c_str3," *", Coeff_Names(3));
hs2_c= strcat(" ",Coeff_tabl_c_str2," PC1",' ',Coeff_tabl_c_str3," PC2"); 
 
hs2_c= strcat(hs2_c,",  r = ",num2str(round(ry1(2)*1000)/1000));
 
 
  Coeff_b_c_str1= num2str(round(b(1,1)*1000)/1000);
  
  hs3_c= strcat(" ",Coeff_b_c_str1," ", Coeff_Names(2));
  hs3_c= strcat(" ",Coeff_b_c_str1," PC1");
  hs3_c= strcat(hs3_c,",  r = ",num2str(round(ry2(2)*1000)/1000));
 
   h_leg=legend([ h1 h2 h3],  hs1_c, hs2_c, hs3_c); 
   
  set(h_leg, 'location', 'NorthOutside','EdgeColor',[0 0 0],'FontSize',16,'FontWeight','bold' ); % ,'color','none'
 
 hl=hline(0);
 set(hl,'LineWidth',3);
 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % letter
  letter='a'; 
  letter_pos=[93 400 50 50];
  
          TextBox = uicontrol('style','text');
          set(TextBox,'String',letter,'position',letter_pos,'FontWeight','bold','FontSize',30 ); % x position y position xsize ysize
          set(TextBox,'foregroundcolor', [0 0 0], ...
         'backgroundcolor', [1 1 1]);
     
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % label    
  letter_pos=[1025 406 200 45];
  
          TextBox = uicontrol('style','text');
          set(TextBox,'String',name_str,'position',letter_pos,'FontWeight','bold','FontSize',30 ); % x position y position xsize ysize
          set(TextBox,'foregroundcolor', [0 0 0], ...
         'backgroundcolor', [1 1 1]);     
     
 %%%%%%%%%%%%%%%%%%%%%%%%%
 
 
 ylabel(['SAT anomaly (',char(176),'C)'],'FontWeight','bold','FontSize',18 );
 
         set(gca,...
          'linewidth',2,...
          'FontWeight','bold' );
 
 
%
%%% Save Fig. 

     filename=['stepwise_MLR','_', name_str]; 
      
      
      filedir ='C:\PHD\matlab_storage_of_output_files\figures\';
      savefilename_c=strcat(filedir,filename);

 % save as png 
orient landscape
% increase quality of picture by increase to -r500


  export_fig('-png','-nocrop','-painters', '-depsc','-opengl', '-r150',savefilename_c); % PNG-nocrop' 
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RICE AWS SAT 

  load('C:\PHD\matlab_storage_of_output_files\RICE_monthly_SAT_anomaly_Ma.mat');
  load('C:\PHD\matlab_storage_of_output_files\ERA-Interim_PCs_2mT_lim0-360_-30_-90_1979-2014_annual.mat');
  name_str='SAT';
            
a_c=385;
b_c=418;

y=MA_SAT_RICE_AWS_monthly((1:34),4);

MA_PCs_save_monthly_c=MA_PCs_save_monthly((a_c:b_c),2);


[X_PC,clim]=annave(MA_PCs_save_monthly_c); % just PC1


[r,p]=corrcoef_df(y,X_PC);

rp=[r(2),p(2)]

p_level(p(2))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%% Stepwise MLR
  
  [X_PC1,clim1]=annave(MA_PCs_save_monthly((a_c:b_c),2));
  [X_PC2,clim2]=annave(MA_PCs_save_monthly((a_c:b_c),3));
  [X_PC3,clim3]=annave(MA_PCs_save_monthly((a_c:b_c),4));
  
  X_PC_step=zeros(34,3);
  X_PC_step(:,1)=X_PC1;
  X_PC_step(:,2)=X_PC2;
  X_PC_step(:,3)=X_PC3;

  LM=stepwiselm(X_PC_step,y,'Upper','linear'); % PCs
  
  
  cur_fitted=LM.Fitted; % Fitted (predicted) response values based on input data

%% Figure RICE AWS 
% Fig. 14b

%plotDiagnostics(LM)

Coeff_tabl=LM.Coefficients;
Coeff_tabl_c=table2array(Coeff_tabl);
Coeff_Names=LM.CoefficientNames;      
      
 X_PC1_cp=X_PC_step(:,1);

b=regress(y,X_PC1_cp);
ry2=corrcoef_df(y,X_PC1_cp*b);


name_str='AWS SAT';

%%%%%%%%%%%%%%%%%%%
% Figure
h=fig('units','inches','width',14,'height',5,'font','Helvetica','fontsize',16,'border','on');

    h1=plot(MA_PCs_save_monthly((a_c:b_c),1),y,'-','LineWidth',3,'MarkerSize',8,'Color',[0.1,0.1,0.1]); 
      hold on 
    
 %%%%%%%%%%%%%%%%%
%    set(gca,'XLim',[2010.6 2014])
     set(gca,'XLim',[2008.9 2015])
%    set(gca,'YLim',[-6 8])
   set(gca,'YLim',[-15 15]) 
      
    line1=cur_fitted;

    h2=plot(MA_PCs_save_monthly((a_c:b_c),1),line1,'--*','LineWidth',3,'MarkerSize',8,'Color',[1,0.3,0.8]);
    
%     h3=plot(MA_PCs_save_monthly((a_c:b_c),1),X_PC1_cp*b,'--','LineWidth',3,'MarkerSize',8,'Color',[0.95,0.45,0.1]);
    
    ry1=corrcoef_df(y,line1);       
% ry1(2)^2
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% legend
 hs1_c='SAT RICE';
 % PCs
 Coeff_tabl_c_str1= num2str(round(Coeff_tabl_c(1,1)*1000)/1000);
 Coeff_tabl_c_str2= num2str(round(Coeff_tabl_c(2,1)*1000)/1000);

 hs2_c= strcat(" ",Coeff_tabl_c_str2," PC1");%,' ',Coeff_tabl_c_str3," *", Coeff_Names(3));
 
 hs2_c= strcat(hs2_c,",  r = ",num2str(round(ry1(2)*1000)/1000));
 
 
 
   h_leg=legend([ h1 h2 ],  hs1_c, hs2_c); 
   
  set(h_leg, 'location', 'NorthOutside','EdgeColor',[0 0 0],'FontSize',16,'FontWeight','bold' ); % ,'color','none'
 
 hl=hline(0);
 set(hl,'LineWidth',3);
 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % letter
  letter='b';
  letter_pos=[93 400 50 50];
  
          TextBox = uicontrol('style','text');
          set(TextBox,'String',letter,'position',letter_pos,'FontWeight','bold','FontSize',30 ); % x position y position xsize ysize
          set(TextBox,'foregroundcolor', [0 0 0], ...
         'backgroundcolor', [1 1 1]);
     
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % label  
  letter_pos=[1025 406 200 45];
  
          TextBox = uicontrol('style','text');
          set(TextBox,'String',name_str,'position',letter_pos,'FontWeight','bold','FontSize',30 ); % x position y position xsize ysize
          set(TextBox,'foregroundcolor', [0 0 0], ...
         'backgroundcolor', [1 1 1]);     
     
 %%%%%%%%%%%%%%%%%%%%%%%%%
 
 
 ylabel(['SAT anomaly (',char(176),'C)'],'FontWeight','bold','FontSize',18 );
 xlabel('Year','FontWeight','bold','FontSize',18 );
 
         set(gca,...
          'linewidth',2,...
          'FontWeight','bold' );

%%% Save Fig. 

     filename=['stepwise_MLR','_', name_str,'_',hs1_c]; 
      
      
      filedir ='C:\PHD\matlab_storage_of_output_files\figures\';
      savefilename_c=strcat(filedir,filename);

% save as png 
orient landscape
% increase quality of picture by increasing, e.g. -r200

export_fig('-png','-nocrop','-painters', '-depsc','-opengl', '-r150',savefilename_c); % PNG-nocrop'     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
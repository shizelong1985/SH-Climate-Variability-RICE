%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stepwise_regression_c.m
% Matlab 2017a
% Daniel Emanuelsson 2018
% subfunctions
% *fig.m fileexchange

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% using the stepwise fitting algorithm.  This method proceeds by adding and removing terms
% in an iterative fitting process, looking at the sum of squared errors (SSE) as a performance
% criterion.   Applying  this  algorithm  allow  us  to  more  confidently  reduce  the  size  of  the  in-
% puts  dataset  by  retaining  only  the  most  influential  ones. 


clear all
close all
%

yr_s=1979;
yr_e=2009;

season='annual';

% Load RICE Annual

% cor_nr=2;  % (1) 1213B (2) winstrup final --use  (3) old
%%%%%%%%%%%%%%%%%

    
      load('C:\PHD\matlab_storage_of_output_files\RICE_combined_Deep_1213B_c19.mat'); % Winstrup 2017 May age-scale
%     load('C:\PHD\matlab_storage_of_output_files\RICE_combined_Deep_1213B_c20.mat'); % Winstrup 2017 Feb age-scale
    stacked_record_annual_Ma=MA_save;
    date_annual=stacked_record_annual_Ma(:,1);


        start_t=find(date_annual==yr_s);
        end_t=find(date_annual==yr_e); 

    y=detrend(stacked_record_annual_Ma((start_t:end_t),2)); % 2 dD 1213B 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Load PCs

param_nr=1; % (1) SIC, c  (2) SAT, b (3) z500  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%

if param_nr==1 % SIC

    load('C:\PHD\matlab_storage_of_output_files\HadISST_PCs_SIC_lim-150--30_-64_-75_1979-2014_annual.mat') 
        

name_str='SIC';
a_c=1; % 1999
b_c=31;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif param_nr==2 % SAT

       load('C:\PHD\matlab_storage_of_output_files\ERA-Interim_PCs_2mT_lim0-360_-30_-90_1979-2014_annual.mat');
       name_str='SAT';
 a_c=1;
%a_c=12; % 1990
%a_c=17; % 1995
b_c=31; 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif param_nr==3 % z500
  
  load('C:\PHD\matlab_storage_of_output_files\ERA-Interim_PCs_z500_lim0-360_-20_-90_1979-2014_annual.mat')
  % MA_PCs_save_z500=MA_PCs_save;
       name_str='z500';
a_c=1;
b_c=31; 




end


 X=[MA_PCs_save((a_c:b_c),(2:4)) ones(length([a_c:b_c]),1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% help stepwiselm
% if param_nr==1 || param_nr==2


%%%%%%%%% stepwiselm %%%%%%%%%%

% Stepwise regression differs from multiple regression in that it starts with just one “best” 
% predictor based on some common goodness of fit criterion such as the sim of squared error (default),
% Akaike information criterion (AIC), or Bayesian information criterion (BIC). Then the algorithm will
% try to add other predictors one at a time until the goodness of fit does not improve. Then it goes 
% in reverse and tries to remove predictors one at a time, again stopping when the goodness of fit does not improve. 

%%%%%%%%%%%%%%%%%
% Matlab function
% LM=stepwiselm(X,y);
% 
% % LM=stepwiselm(X,y,'Upper');
% LM=stepwiselm(X,y,'linear');
% LM=stepwiselm(X,y,'linear','Upper','linear');
% LM=stepwiselm(X,y,'linear','Lower','linear');


method_nr=1; % (1) default, (2)- relaxed p-value

if method_nr==1
    LM=stepwiselm(X,y,'Upper','linear'); % PCs
elseif method_nr==2
    LM=stepwiselm(X,y,'Upper','linear','PEnter',0.08); % PCs
end



% The input parameter name/value pair "Upper" with its value "linear"
% indicates that we only want linear fucntions of the indpenedent variables
% as predictors, not products of them.


b_Coeff=LM.CoefficientCovariance;

cur_fitted=LM.Fitted; % Fitted (predicted) response values based on input data

% Which constants does it use

% plotDiagnostics(LM)

Coeff_tabl=LM.Coefficients;
Coeff_tabl_c=table2array(Coeff_tabl);
Coeff_Names=LM.CoefficientNames;
%     if param_nr==2
%     Coeff_tabl1=LM1.Coefficients;
%     Coeff_tabl_c1=table2array(Coeff_tabl1);
%     Coeff_Names1=LM1.CoefficientNames;
%     
%     cur_fitted_c1=LM1.Fitted;
%     end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure

h=fig('units','inches','width',14,'height',5,'font','Helvetica','fontsize',16,'border','on');

    h1=plot([yr_s:2009],y,'-','LineWidth',3,'MarkerSize',8,'Color',[0.1,0.1,0.1]); 
      hold on  
    
  if param_nr==2 % SAT
   letter='b';
  elseif param_nr==1 %SIC
   letter='c';
  elseif param_nr==3 %z500
   letter='a'; 
  elseif param_nr==4 %SST
   letter='a'; 
  end
  
  %%%%%%%%%%%%%%%%%
  set(gca,'XLim',[yr_s-1 2010])
  set(gca,'YLim',[-30 35]) 

 if param_nr==2
    

alt_nr=2; % (1/2) should be the same

    if alt_nr==1
        line1=Coeff_tabl_c(1,1)+Coeff_tabl_c(2,1)*X(:,1)+Coeff_tabl_c(3,1)*X(:,8);
%         line_c1=Coeff_tabl_c1(1,1)+Coeff_tabl_c1(2,1)*X(:,1); 
    elseif alt_nr==2
        
        line1=cur_fitted;
%         line_c1=cur_fitted_c1;
    end 
 

% h2=plot([1979:2009],line_c1,'--*','LineWidth',3,'MarkerSize',8,'Color',[0.8,0.3,0.2]);
 h2=plot([1979:2009],line1,'--*','LineWidth',3,'MarkerSize',8,'Color',[0.65,0.05,0.0]); 



% PC1-10 
ry1=corrcoef_df(y,line1);
% PC1-4     
% ry_c1=corrcoef_df(y,line_c1);
% ry1(2)^2

%%%%%%%%%%%%   
% legend
  hs1_c='dD detrended';
 % PCs
 Coeff_tabl_c_str1= num2str(round(Coeff_tabl_c(1,1)*1000)/1000);
 Coeff_tabl_c_str2= num2str(round(Coeff_tabl_c(2,1)*1000)/1000);
%  Coeff_tabl_c_str3= num2str(round(Coeff_tabl_c(3,1)*1000)/1000);
  
  
%  hs2_c= strcat(" ",Coeff_tabl_c_str1," ",Coeff_tabl_c_str2," ", Coeff_Names(2));%,' +',Coeff_tabl_c_str3," *", Coeff_Names(3));
 hs2_c= strcat(" ",Coeff_tabl_c_str1," ",Coeff_tabl_c_str2," ", 'PC1');%
 hs2_c= strcat(hs2_c,",  r = ",num2str(round(ry1(2)*1000)/1000));
   
%      h_leg=legend([ h1 h2 h3],  hs1_c, hs2_c, hs3_c); 
  h_leg=legend([ h1 h2],  hs1_c, hs2_c); 
   
  set(h_leg, 'location', 'NorthOutside','EdgeColor',[0 0 0],'FontSize',18,'FontWeight','bold' ); % ,'color','none'
 
 hl=hline(0);
 set(hl,'LineWidth',3);
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 % SIC
elseif param_nr==1
 
alt_nr=2; % (1/2) should be the same

    if alt_nr==1
        line1=Coeff_tabl_c(1,1)+Coeff_tabl_c(2,1)*X(:,1)+Coeff_tabl_c(3,1)*X(:,2);
    elseif alt_nr==2
        line1=cur_fitted;
    end 


h2=plot([1979:2009],line1,'--*','LineWidth',3,'MarkerSize',8,'Color',[0.6,0.0,0.4]);    
    
 % PC
ry1=corrcoef_df(y,line1);
 
%%%%%%%%%%%%   
% legend

 % PC 
 Coeff_tabl_c_str1= num2str(round(Coeff_tabl_c(1,1)*1000)/1000);
 Coeff_tabl_c_str2= num2str(round(Coeff_tabl_c(2,1)*1000)/1000);
 Coeff_tabl_c_str3= num2str(round(Coeff_tabl_c(3,1)*1000)/1000);
  
 hs1_c='dD detrended'; 
%  hs2_c= strcat(" ",Coeff_tabl_c_str1," ",Coeff_tabl_c_str3," ", Coeff_Names(3)," ",Coeff_tabl_c_str2," ", Coeff_Names(2)); 
 hs2_c= strcat(" ",Coeff_tabl_c_str1," ",Coeff_tabl_c_str3," ", 'PC2'," ",Coeff_tabl_c_str2," ", 'PC1');  
 
 hs2_c= strcat(hs2_c,",  r = ",num2str(round(ry1(2)*1000)/1000));
 
      h_leg=legend([ h1 h2 ],  hs1_c, hs2_c); 

   
  set(h_leg, 'location', 'NorthOutside','EdgeColor',[0 0 0],'FontSize',18,'FontWeight','bold' ); % ,'color','none'
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  elseif param_nr==3 % z500
      
        alt_nr=2; % (1/2) should be the same

    if alt_nr==1
        line1=Coeff_tabl_c(1,1)+Coeff_tabl_c(2,1)*X(:,1)+Coeff_tabl_c(3,1)*X(:,8);
    elseif alt_nr==2
        line1=cur_fitted;
    end      
      
    h2=plot([1979:2009],line1,'--*','LineWidth',3,'MarkerSize',8,'Color',[0.0,0.2,0.6]); 
    
      
   % PCs 
    ry1=corrcoef_df(y,line1);
 
 %%%%%%%%%%%%   
% legend
    Coeff_tabl_c_str1= num2str(round(Coeff_tabl_c(1,1)*1000)/1000);
    Coeff_tabl_c_str2= num2str(round(Coeff_tabl_c(2,1)*1000)/1000);
    Coeff_tabl_c_str3= num2str(round(Coeff_tabl_c(3,1)*1000)/1000);
  
    hs1_c='dD detrended'; 
    
%     hs2_c= strcat(" ",Coeff_tabl_c_str1," +",Coeff_tabl_c_str3," ", Coeff_Names(3),...
%         " ",Coeff_tabl_c_str2," ", Coeff_Names(2));
    hs2_c= strcat(" ",Coeff_tabl_c_str1," +",Coeff_tabl_c_str3," ", 'PC2',...  % PC instead of x
        " ",Coeff_tabl_c_str2," ", 'PC1');
    
    hs2_c= strcat(hs2_c,",  r = ",num2str(round(ry1(2)*1000)/1000));
 
    h_leg=legend([ h1 h2 ],  hs1_c, hs2_c); 

   
    set(h_leg, 'location', 'NorthOutside','EdgeColor',[0 0 0],'FontSize',18,'FontWeight','bold' ); % ,'color','none'    
  
      
 end

 
 hl=hline(0);
 set(hl,'LineWidth',3);
 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Letter
  letter_pos=[93 400 50 50];
  
          TextBox = uicontrol('style','text');
          set(TextBox,'String',letter,'position',letter_pos,'FontWeight','bold','FontSize',30 ); % x position y position xsize ysize
          set(TextBox,'foregroundcolor', [0 0 0], ...
         'backgroundcolor', [1 1 1]);
     
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % label
 % parameter string 
  letter_pos=[1125 400 100 45];
  
          TextBox = uicontrol('style','text');
          set(TextBox,'String',name_str,'position',letter_pos,'FontWeight','bold','FontSize',30 ); % x position y position xsize ysize
          set(TextBox,'foregroundcolor', [0 0 0], ...
         'backgroundcolor', [1 1 1]);     
     
 %%%%%%%%%%%%%%%%%%%%%%%%%
 ylabel(['{\delta}D anomaly (',char(8240),')'],'FontWeight','bold','FontSize',18 );
 
         set(gca,...
          'linewidth',2,...
          'FontWeight','bold' );
          
  if param_nr==1
         xlabel('Year','FontWeight','bold','FontSize',18 );
  end
 
%
%%% Save Fig. 
     filename=['stepwise_MLR','_', name_str];  
      filedir ='C:\PHD\matlab_storage_of_output_files\figures\';
      savefilename_c=strcat(filedir,filename);

 % save as png 
orient landscape
% increase quality of picture by increase to -r300
 export_fig('-png','-nocrop','-painters', '-depsc','-opengl', '-r170',savefilename_c); % PNG-nocrop' 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%% SIC ADP
  
       load('C:\PHD\matlab_storage_of_output_files\Ma_ADP_index_HadISST_SIC.mat');
       name_str='SIC ADP';
       
a_c=8;
b_c=38;  

    
X=[MA_ADP_save((a_c:b_c),(2)) ones(31,1)];
    
b=regress(y,X);

LM=stepwiselm(X,y,'Upper','linear'); % PCs


Coeff_tabl=LM.Coefficients;
Coeff_tabl_c=table2array(Coeff_tabl);
  
b_Coeff=LM.CoefficientCovariance;

cur_fitted=LM.Fitted; % Fitted (predicted) response values based on input data


Coeff_tabl=LM.Coefficients;
Coeff_tabl_c=table2array(Coeff_tabl);
Coeff_Names=LM.CoefficientNames;  
  
 line1=cur_fitted;
 
 
    % PCs 
    ry1=corrcoef_df(y,line1);

    
    %ry1(2)^2;
%

h=fig('units','inches','width',14,'height',5,'font','Helvetica','fontsize',16,'border','on');

    h1=plot([1979:2009],y,'-','LineWidth',3,'MarkerSize',8,'Color',[0.1,0.1,0.1]); 
  

    hold on

    h2=plot([1979:2009],line1,'-x','LineWidth',3,'MarkerSize',8,'Color',[0.62,0.52,0.52]); % SIC ADP
    
%       axestext_c(1.15,1.00,name_str,'FontWeight','bold','FontSize',20 );   

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % letter
  letter_pos=[90 400 50 50];
  letter='b';     
 
          TextBox = uicontrol('style','text');
          set(TextBox,'String',letter,'position',letter_pos,'FontWeight','bold','FontSize',30 ); % x position y position xsize ysize
          set(TextBox,'foregroundcolor', [0 0 0], ...
         'backgroundcolor', [1 1 1]);
%%%%%%%%%%%%%%%%%%%%%%%%%
% label     
  letter_pos=[1005 400 250 45];
  
          TextBox = uicontrol('style','text');
          set(TextBox,'String',name_str,'position',letter_pos,'FontWeight','bold','FontSize',30 ); % x position y position xsize ysize
          set(TextBox,'foregroundcolor', [0 0 0], ...
         'backgroundcolor', [1 1 1]);     
     
 %%%%%%%%%%%%%%%%%%%%%%%%%

         set(gca,'XLim',[1978 2010])
         set(gca,'YLim',[-30 35])
         %

 %%%%%%%%%%%%%%%%%%%%%%%%%
 % legend
    hs2='ADP index';   
    hs1_c='dD detrended';
%     hs2_c=[hs2,', ','r = ',num2str(0.001*round(ry1(1,2)*1000))];



    Coeff_tabl_c_str1= num2str(round(Coeff_tabl_c(1,1)*1000)/1000);
    Coeff_tabl_c_str2= num2str(round(Coeff_tabl_c(2,1)*1000)/1000);


%      hs2_c= strcat(" ",Coeff_tabl_c_str1," ",Coeff_tabl_c_str2," ", Coeff_Names(2));%,' +',Coeff_tabl_c_str3," *", Coeff_Names(3));
     hs2_c= strcat(" ",Coeff_tabl_c_str1," ",Coeff_tabl_c_str2," ", 'ADP');
     hs2_c= strcat(hs2_c,",  r = ",num2str(round(ry1(2)*1000)/1000));  
     h_leg=legend([ h1 h2],  hs1_c, hs2_c); 

   
  set(h_leg, 'location', 'NorthOutside','EdgeColor',[0 0 0],'FontSize',18 ,'FontWeight','bold'); % ,'color','none'

       xlabel('Year','FontWeight','bold','FontSize',18 );
       ylabel(['{\delta}D anomaly (',char(8240),')'],'FontWeight','bold','FontSize',18 );
       
        hl=hline(0);
        set(hl,'LineWidth',3);

         set(gca,...
          'linewidth',2,...
          'FontWeight','bold' );

%%% Save Fig. 

     filename=['MLR','_', name_str]; 
      
      
      filedir ='C:\PHD\matlab_storage_of_output_files\figures\';
      savefilename_c=strcat(filedir,filename);

 % save as png 
orient landscape
% increase quality of picture by increase to -r500
export_fig('-png','-nocrop','-painters', '-depsc','-opengl', '-r170',savefilename_c); % PNG-nocrop' 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  %% All PCs
  %clear all
  
  %load('C:\PHD\matlab_storage_of_output_files\HadISST_PCs_SIC_lim-150--30_-65_-90_1979-2014_annual.mat')
  load('C:\PHD\matlab_storage_of_output_files\HadISST_PCs_SIC_lim-150--30_-64_-75_1979-2014_annual.mat')
  MA_PCs_save_SIC=MA_PCs_save;

  load('C:\PHD\matlab_storage_of_output_files\ERA-Interim_PCs_2mT_lim0-360_-30_-90_1979-2014_annual.mat');
  MA_PCs_save_SAT=MA_PCs_save;

  load('C:\PHD\matlab_storage_of_output_files\ERA-Interim_PCs_z500_lim0-360_-20_-90_1979-2014_annual.mat')
  MA_PCs_save_z500=MA_PCs_save;
  
  
%%

alt_nr=1;%%%%%% start year (1) 1979 use

if alt_nr==1
    a_c=1;
elseif alt_nr==2
    a_c=23;
elseif alt_nr==3
    a_c=13;
end

b_c=31; 
% doesnt matter which order you put them in 

PC_feed_nr=1; % (1) All PCs (2) two first predictors

if PC_feed_nr==1
  X=[ MA_PCs_save_SAT((a_c:b_c),(2:4)) MA_PCs_save_SIC((a_c:b_c),(2:4))  MA_PCs_save_z500((a_c:b_c),(2:4))  ones(length([a_c: b_c]),1)];
elseif PC_feed_nr==2 
  X=[ MA_PCs_save_SAT((a_c:b_c),(2)) MA_PCs_save_z500((a_c:b_c),(3))  ones(31,1)];
end
  
  
LM=stepwiselm(X,y(a_c:b_c,1),'Upper','linear'); % All PCs

  
b_Coeff=LM.CoefficientCovariance;

cur_fitted=LM.Fitted; % Fitted (predicted) response values based on input data

% Which constants does it use

% plotDiagnostics(LM)

Coeff_tabl=LM.Coefficients;
Coeff_tabl_c=table2array(Coeff_tabl);
Coeff_Names=LM.CoefficientNames;


      line1=cur_fitted;
 
 
    % PCs 
    ry1=corrcoef_df(y(a_c:b_c,1),line1);
 
    %ry1(2)^2;

    Coeff_tabl_c_str_int= num2str(round(Coeff_tabl_c(1,1)*1000)/1000); % intercept
    Coeff_tabl_c_str1= num2str(round(Coeff_tabl_c(2,1)*1000)/1000); % x1
    Coeff_tabl_c_str2= num2str(round(Coeff_tabl_c(3,1)*1000)/1000); % x2
    Coeff_tabl_c_str3= num2str(round(Coeff_tabl_c(4,1)*1000)/1000); % x3
    Coeff_tabl_c_str8= num2str(round(Coeff_tabl_c(5,1)*1000)/1000); % x8
    


    hs2_c= strcat(Coeff_tabl_c_str_int," ",Coeff_tabl_c_str1," CP1SAT +",Coeff_tabl_c_str8," CP2z500",...
         " +",Coeff_tabl_c_str2," CP2SAT", Coeff_tabl_c_str3," CP3SAT" );
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%% Figure
    h=fig('units','inches','width',14,'height',5,'font','Helvetica','fontsize',16,'border','on');

    h1=plot([1979:2009],y,'-','LineWidth',3,'MarkerSize',8,'Color',[0.1,0.1,0.1]);
    hold on

    h2=plot([1979:2009],line1,'-x','LineWidth',3,'MarkerSize',8,'Color',[0.8,0.1,0.32]); 
 
  %%%%%%%%%%%%%
  % legend
       hs1_c='dD detrended';
       hs2_c= strcat(hs2_c,",  r = ",num2str(round(ry1(2)*1000)/1000));
 
    h_leg=legend([ h1 h2 ],  hs1_c, hs2_c); 

   
    set(h_leg, 'location', 'NorthOutside','EdgeColor',[0 0 0],'FontSize',18,'FontWeight','bold' ); % ,'color','none'    
 %%%%%%%%%%%%
 
           set(gca,'XLim',[1978 2010])
%          set(gca,'YLim',[-30 35])

 name_str='All PCs';    
 hl=hline(0);
 set(hl,'LineWidth',3);
 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % letter
  letter_pos=[93 400 50 50];
  letter='a';
          TextBox = uicontrol('style','text');
          set(TextBox,'String',letter,'position',letter_pos,'FontWeight','bold','FontSize',30 ); % x position y position xsize ysize
          set(TextBox,'foregroundcolor', [0 0 0], ...
         'backgroundcolor', [1 1 1]);
     
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % label    
  letter_pos=[1025 310 190 45];
  
          TextBox = uicontrol('style','text');
          set(TextBox,'String',name_str,'position',letter_pos,'FontWeight','bold','FontSize',30 ); % x position y position xsize ysize
          set(TextBox,'foregroundcolor', [0 0 0], ...
         'backgroundcolor', [1 1 1]);     
     
 %%%%%%%%%%%%%%%%%%%%%%%%%
 ylabel(['{\delta}D anomaly (',char(8240),')'],'FontWeight','bold','FontSize',18 );
%  xlabel('Year','FontWeight','bold','FontSize',18 ); 
         set(gca,...
          'linewidth',2,...
          'FontWeight','bold' );
      
%%% Save Fig.
     filename=['stepwise_MLR','_', name_str]; 
      
      
      filedir ='C:\PHD\matlab_storage_of_output_files\figures\';
      savefilename_c=strcat(filedir,filename);

 % save as png 
orient landscape
% increase quality of picture by increase to -r300

  export_fig('-png','-nocrop','-painters', '-depsc','-opengl', '-r170',savefilename_c); % PNG-nocrop' 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
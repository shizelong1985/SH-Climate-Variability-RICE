%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RICE_d18O_anomaly_EOF_paper
% d18O annual time series
% DE 2020
% MATLAB 2018a
%%%%%%%%%%%%%%%%%%%%%%%%%%
% uses imput from composite_SAM_PSA1_in_phase_ERA_I_z500_2mT_HadISST_SIC_c15.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

load('C:\PHD\matlab_storage_of_output_files\RICE_combined_Deep_1213B_c23.mat');

f1=fig('units','inches','width',18,'height',5,'font','Helvetica','fontsize',18,'border','on'); 
% set(gca,'YLim',[-3 4])
    
      hold on
      box on
      
      
%plot(MA_save(:,1),anomaly(nanmoving_average(MA_save(:,3),10)),'-k','LineWidth',2.5)
plot(MA_save((1980:end),1),MA_save((1980:end),3),'-k','LineWidth',2.5)

plot(MA_save((1980:end),1),MA_save((1980:end),3),'.k','MarkerSize',22)

wr_c=nanmean(MA_save((1980:end),3));
lrc3=hline(wr_c,':');
  
set(lrc3,'Color',[0.1,0.1,0.1],'LineWidth',1.5);
  

wr_c1=nanmean(MA_save((1980:end),3))+nanstd(MA_save((1980:end),3)); % std lines
lrc4=hline(wr_c1,':');
set(lrc4,'Color',[0.9,0.0,0.0],'LineWidth',1.5);  


wr_c1=nanmean(MA_save((1980:end),3))-nanstd(MA_save((1980:end),3)); % std lines
lrc4=hline(wr_c1,':');
set(lrc4,'Color',[0.9,0.0,0.0],'LineWidth',1.5); 


yr_s=1979;
yr_e=2011;

set(gca,'XLim',[1978 2012])
set(gca,'YLim',[-27 -19])
%set(gca,'XLim',[1900 2011])

t_c=[1979:2011]';

% from composite_SAM... file

ind_c1 =[1
    7
    11
    20
    21
    30
    32];


plot(t_c(ind_c1,1),wr_c,'*k','MarkerSize',14 ) % years



% % rotate
ind_c2 =[ 2
    13
    24
    31];



plot(t_c(ind_c2,1),wr_c,'ok','MarkerSize',14 ) % years

 
              set(gca,...
          'linewidth',2,...
          'FontWeight','bold' );

     xlabel('Year','FontWeight','bold','FontSize',20 );
     ylabel(['{\delta}^1^8O (',char(8240),')'],'FontWeight','bold','FontSize',20 );
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [x y width height] for the figure, increase margins
    pos1 = get(gca, 'Position');
    pos1(1) = 0.2;
    pos1(3) = 0.6;
    set(gca, 'Position', pos1)
    
%%%%%%%%%%%%
 hs1=['{\delta}^1^8O']; %(',char(8240),')'];   

l2=legend(['RI ', hs1]);
 set(l2,'FontSize',14, 'FontWeight', 'bold', 'location', 'NorthWest');    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% letter

letter='a';

TextBox = uicontrol('style','text');
letter_position=   [240 380 60 60];
   
set(TextBox,'String',letter,'position',letter_position,'FontWeight','bold','FontSize',40 ); 
                 set(TextBox,'foregroundcolor', [0 0 0], ...
            'backgroundcolor', [1 1 1]);
%%%%%%%%%%%%    
    
     
     
%%%%%%%%%%%%%%%%          
save_fig=1;

if save_fig==1
iso_str='d18O';

filename=['RI_anomaly_',iso_str,'_',num2str(yr_s),'-',num2str(yr_e),'_fig'];
    
     %edit  save_fig_DE
filedir ='C:\PHD\matlab_storage_of_output_files\figures\';
savefilename_c=strcat(filedir,filename);

 % save as png 
orient landscape
% increase quality of picture by increase to -r500


 export_fig('-png','-nocrop','-painters', '-depsc','-opengl', '-r190',savefilename_c);            
 export_fig('-pdf','-nocrop','-painters', '-depsc','-opengl', '-r190',savefilename_c);          
 
% export_fig(savefilename_c,'-transparent','-eps','-painters');
 
 
end            
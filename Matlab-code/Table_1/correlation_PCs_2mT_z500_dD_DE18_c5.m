%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% correlation_PCs_2mT_z500_dD_DE18_c5.m
%%%% For Table 1
% Daniel Emanuelsson 2018
% Matlbab 2017a
% Github version 1  11-09-2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subfunctions
% * corrcoef_df.m   UoW Steig
% * p_level.m       DE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iso=0; % (1/0) dD on/off
PC_pos=2+1; %(1-PC1, 2-PC2, 3-PC3)%%

PCs_nr=1; %(1 z500 / 2 SAT / 3 SIC) =-
PCs_nr2=1;%(1 SAT / 2 SIC)

%%%%%%%%%%% for z500 SAT corr corresponding pattern pairs . Comment
%%%%%%%%%%% uncomment to get parantheses calc.
if PCs_nr==1 && PCs_nr2==1 && iso==0
    if PC_pos==3
        PC_pos_c=PC_pos+1;
    elseif PC_pos==4 
        PC_pos_c=PC_pos-1;
    elseif PC_pos==2
        PC_pos_c=PC_pos;
    end
else
   PC_pos_c=PC_pos;   
end
 
%%%%%%%%%%%%%%%%%%%%%%%
if PCs_nr==1
        PCs_st='z500';
        lat1=-20;
elseif PCs_nr==2
        PCs_st='2mT';
        lat1=-30;
elseif PCs_nr==3
        PCs_st='SIC';    
end

    if iso==1
    
        yr_c=31;
    else
        yr_c=35;
    end

    %(1)
      
       if (PCs_nr==1 || PCs_nr==2)
            load(['C:\PHD\matlab_storage_of_output_files\ERA-Interim_PCs_',PCs_st,'_lim0-360_',num2str(lat1),'_-90_1979-2014_annual.mat']);    
       elseif PCs_nr==3           
            load('C:\PHD\matlab_storage_of_output_files\HadISST_PCs_SIC_lim-150--30_-64_-75_1979-2014_annual.mat'); % SIC
       end
        
    
    % Same as fig
if PCs_nr==1 % z500  
    if PC_pos==2 % SAM
        
        factor_nr=1;
        
    elseif PC_pos==3 % PSA1  [polarity of PSA patterns defined as in Kidson 1988, Fig 4b,c]
        
         factor_nr=1;
         
    elseif PC_pos==4 % PSA2
        
         factor_nr=-1;
    end
    
elseif PCs_nr==2 % SAT 2mT
    
    % Same as regres fig
    if PC_pos==2 % PC1
         factor_nr=1;
    elseif PC_pos==3 % PC2
         factor_nr=1;
    elseif PC_pos==4 % PC3
        factor_nr=1;
    end
        
  elseif PCs_nr==3 % SIC   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    
          % Same as regres fig
    if PC_pos==2 % PC1
         factor_nr=1;
    elseif PC_pos==3 % PC2
         factor_nr=-1;
    elseif PC_pos==4 % PC3
        factor_nr=-1;   
    end  
end
    
X=detrend(MA_PCs_save((1:yr_c),PC_pos))* factor_nr;


    if iso==0
          
           if PCs_nr2==1 %SAT
              load('C:\PHD\matlab_storage_of_output_files\ERA-Interim_PCs_2mT_lim0-360_-30_-90_1979-2014_annual.mat'); 
              
                if PC_pos==2 % PC1 SAT
                factor_nr= 1;
                elseif PC_pos==3 % PC2
                factor_nr= 1;
                elseif PC_pos==4 % PC3
                factor_nr= 1;
                end
              
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
           elseif PCs_nr2==2
            
             load('C:\PHD\matlab_storage_of_output_files\HadISST_PCs_SIC_lim-150--30_-64_-75_1979-2014_annual.mat'); % SIC
                if PC_pos==2 % PC1
                    factor_nr=1;
                elseif PC_pos==3 % PC2
                    factor_nr=-1;
                elseif PC_pos==4 % PC3
                    factor_nr=-1;
                end                  
           end

    Y=detrend(MA_PCs_save((1:yr_c),PC_pos_c))* factor_nr;
        
    elseif iso==1

        load('C:\PHD\matlab_storage_of_output_files\RICE_combined_Deep_1213B_c19.mat'); % Winstrup 2017 May age-scale
        date_annual=MA_save(:,1);
        stacked_record_annual_Ma=MA_save;
        Y=detrend(stacked_record_annual_Ma((80:110),2));
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%correlation
%  factor_nr
[r,p]=corrcoef_df(X,Y);
rp=[r(2),p(2)]
p_level(p(2))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
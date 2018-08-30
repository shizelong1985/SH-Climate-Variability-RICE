function IPO_vline(wr)
%%%%%%%%%%%%%%%%%%%
% IPO vertical plot lines and labels
% subfunction
% *vline.m          fileexchange
% *axestext_c.m     JR code
%%%%%%%%%%%%%%%%%%

    if wr==7 || wr==8
    
    x=[1978 2001]; % as in England et al. 2014
    
   elseif wr==9 || wr==10 || wr==15 
       x=[1979    
       2000];     
      

   elseif wr==11
   
         x=[1659
             1681
             1696
             1713
             1727
             1738
             1746
             1786
             1800
             1823
             1853
             1899
             1909
             1922
             1941      
             1979    
             2000];
       
  elseif wr==12     
   
       x=[1978    
       1999];
   
   elseif wr==14  
   
            x=[1659 % same as 11 but no labels
             1681
             1696
             1713
             1727
             1738
             1746
             1786
             1800
             1823
             1853
             1899
             1909
             1922
             1941      
             1979    
             2000];
    
    else
    x=[1910 1925 1944 1958 1968 1978 1999 2006 2013 ]';

    end


b=size(x);



 lr1=vline(x,'--');% IPO
 set(lr1,'Color',[0.2,0.9,0.2],'LineWidth',3.5);
 
         if wr==1
            axestext_c(0.94+0.05 ,0.005,'IPO-','FontWeight','bold','FontSize',20 );
            axestext_c(0.83+0.03 ,0.005,'IPO+','FontWeight','bold','FontSize',20 );
            axestext_c(0.63-0.08 ,0.005,'IPO-','FontWeight','bold','FontSize',20 );
            axestext_c(0.46-0.18 ,0.005,'IPO+','FontWeight','bold','FontSize',20 );
            axestext_c(0.37-0.30 ,0.005,'neu.','FontWeight','bold','FontSize',20 );
            axestext_c(0.30-0.31 ,0.005,'IPO+','FontWeight','bold','FontSize',20 );
            
         elseif wr==2
             adj_nr=0.01;
            axestext_c(0.94+0.023 ,0.005+adj_nr,'IPO-','FontWeight','bold','FontSize',20 );
            axestext_c(0.92+0.00 ,0.005+adj_nr,'neu.','FontWeight','bold','FontSize',20 );
            axestext_c(0.83+0.00 ,0.005+adj_nr,'IPO+','FontWeight','bold','FontSize',20 );
            axestext_c(0.63-0.06 ,0.005+adj_nr,'IPO-','FontWeight','bold','FontSize',20 );
            axestext_c(0.63+0.01 ,0.005+adj_nr,'neu.','FontWeight','bold','FontSize',20 ); 
            axestext_c(0.63+0.07 ,0.005+adj_nr,'IPO-','FontWeight','bold','FontSize',20 );         
            axestext_c(0.46+0.01 ,0.005+adj_nr,'IPO+','FontWeight','bold','FontSize',20 );
            axestext_c(0.36+0.006 ,0.005+adj_nr,'IPO-/neu.','FontWeight','bold','FontSize',20);
            axestext_c(0.30-0.02 ,0.005+adj_nr,'IPO+','FontWeight','bold','FontSize',20 );     % 1900  
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         elseif wr==3
             adj_nr=-0.005;   
            axestext_c(0.94+0.01 ,0.005+adj_nr,'IPO-','FontWeight','bold','FontSize',20 );
            axestext_c(0.92-0.03 ,0.005+adj_nr,'neu.','FontWeight','bold','FontSize',20 );
            axestext_c(0.83-0.05 ,0.005+adj_nr,'IPO+','FontWeight','bold','FontSize',20 );          
            axestext_c(0.63-0.20 ,0.005+adj_nr,'IPO-','FontWeight','bold','FontSize',20 );       
            axestext_c(0.63-0.10 ,0.005+adj_nr,'neu.','FontWeight','bold','FontSize',20 );
            axestext_c(0.63-0.03 ,0.005+adj_nr,'IPO-','FontWeight','bold','FontSize',20 );  
            axestext_c(0.46-0.18 ,0.005+adj_nr,'IPO+','FontWeight','bold','FontSize',20 );
            axestext_c(0.36-0.27 ,0.005+adj_nr,'IPO-/neu.','FontWeight','bold','FontSize',20 ); % 244
            axestext_c(0.30-0.30 ,0.005+adj_nr,'IPO+','FontWeight','bold','FontSize',20 );   % 28
            
         elseif wr==4  % ITASE
             adj_nr=0.01;
             adj_nr2=0.21;
            axestext_c(0.94+0.03 ,0.008+adj_nr2,'IPO-','FontWeight','bold','FontSize',20 );
            axestext_c(0.92+0.00 ,0.005+adj_nr2,'neu.','FontWeight','bold','FontSize',20 );
            axestext_c(0.83+0.00 ,0.005+adj_nr2,'IPO+','FontWeight','bold','FontSize',20 );
            axestext_c(0.63+0.07 ,0.005+adj_nr,'IPO-','FontWeight','bold','FontSize',20 );
            axestext_c(0.63+0.01 ,0.005+adj_nr,'neu.','FontWeight','bold','FontSize',20 );
            axestext_c(0.63-0.06 ,0.005+adj_nr,'IPO-','FontWeight','bold','FontSize',20 );
            axestext_c(0.46+0.01 ,0.005+adj_nr,'IPO+','FontWeight','bold','FontSize',20 );
            axestext_c(0.36+0.00 ,0.005+adj_nr,'IPO-/neu.','FontWeight','bold','FontSize',20 );
            axestext_c(0.30-0.02 ,0.005+adj_nr,'IPO+','FontWeight','bold','FontSize',20 );  % 1900
            
           elseif wr==5   % z500 AS
                            adj_nr=0.01;
                            adj_nr2=0.21;
               
%              axestext_c(0.94+0.06 ,0.005+adj_nr,'IPO-','FontWeight','bold','FontSize',20 );
% add later in combine code
            axestext_c(0.92+0.045 ,0.005+adj_nr,'neu.','FontWeight','bold','FontSize',20 );
            axestext_c(0.83+0.01 ,0.005+adj_nr,'IPO+','FontWeight','bold','FontSize',20 );
            axestext_c(0.63+0.03 ,0.005+adj_nr,'IPO-','FontWeight','bold','FontSize',20 );
            axestext_c(0.63-0.05 ,0.005+adj_nr,'neu.','FontWeight','bold','FontSize',20 );
            axestext_c(0.63-0.17 ,0.005+adj_nr,'IPO-','FontWeight','bold','FontSize',20 );
            axestext_c(0.46-0.15 ,0.005+adj_nr,'IPO+','FontWeight','bold','FontSize',20 );
            axestext_c(0.36-0.27 ,0.005+adj_nr,'IPO-/neu.','FontWeight','bold','FontSize',20 ); % 244
            axestext_c(0.30-0.30 ,0.005+adj_nr,'IPO+','FontWeight','bold','FontSize',20 );   % 28
            
            elseif wr==6              
                    adj_nr=-0.005;
            axestext_c(0.94-0.01 ,0.005+adj_nr,'IPO-','FontWeight','bold','FontSize',20 );
            axestext_c(0.92-0.05 ,0.005+adj_nr,'neu.','FontWeight','bold','FontSize',20 );
            axestext_c(0.83-0.05 ,0.005+adj_nr,'IPO+','FontWeight','bold','FontSize',20 ); 
            
             elseif wr==8             
                    adj_nr=-0.005;
            axestext_c(0.999+0.01 ,0.005+adj_nr,'IPO-','FontWeight','bold','FontSize',20 );
            axestext_c(0.60 ,0.005+adj_nr,'IPO+','FontWeight','bold','FontSize',20 ); 
            
            
             elseif wr==9 % || wr==12            
                    adj_nr=-0.005;
                    
            axestext_c(0.94 ,0.005+adj_nr,'IPO-','FontWeight','bold','FontSize',20 );
            axestext_c(0.60 ,0.005+adj_nr,'IPO+','FontWeight','bold','FontSize',20 );  
            axestext_c(0.15 ,0.005+adj_nr,'IPO-','FontWeight','bold','FontSize',20 );          
            
            elseif wr==11            
            adj_nr=-0.005;
                    
            
%             axestext_c(0.92-0.05 ,0.005+adj_nr,'neu.','FontWeight','bold','FontSize',20 );
            axestext_c(-0 ,0.005+adj_nr,'+','FontWeight','bold','FontSize',20 );
            axestext_c(0.03 ,0.005+adj_nr,'-','FontWeight','bold','FontSize',20 );
            axestext_c(0.08 ,0.005+adj_nr,'+','FontWeight','bold','FontSize',20 );
            axestext_c(0.12 ,0.005+adj_nr,'-','FontWeight','bold','FontSize',20 );
            axestext_c(0.16 ,0.005+adj_nr,'+','FontWeight','bold','FontSize',20 );
            axestext_c(0.19 ,0.005+adj_nr,'-','FontWeight','bold','FontSize',20 );
            axestext_c(0.215 ,0.005+adj_nr,'+','FontWeight','bold','FontSize',20 );
            axestext_c(0.30 ,0.005+adj_nr,'-','FontWeight','bold','FontSize',20 );
            axestext_c(0.36 ,0.005+adj_nr,'+','FontWeight','bold','FontSize',20 );
            axestext_c(0.40 ,0.005+adj_nr,'-','FontWeight','bold','FontSize',20 );
            axestext_c(0.47 ,0.005+adj_nr,'+','FontWeight','bold','FontSize',20 );
            axestext_c(0.55 ,0.005+adj_nr,'-','FontWeight','bold','FontSize',20 );
            axestext_c(0.635 ,0.005+adj_nr,'+','FontWeight','bold','FontSize',20 );
            axestext_c(0.66 ,0.005+adj_nr,'-','FontWeight','bold','FontSize',20 ); 
            axestext_c(0.70 ,0.005+adj_nr,'+','FontWeight','bold','FontSize',20 ); 
            axestext_c(0.78 ,0.005+adj_nr,'-','FontWeight','bold','FontSize',20 ); 
            axestext_c(0.87 ,0.005+adj_nr,'+','FontWeight','bold','FontSize',20 ); 
            axestext_c(0.90 ,0.005+adj_nr,'-','FontWeight','bold','FontSize',20 );
            elseif wr==15             
                    adj_nr=-0.005;
            axestext_c(0.999+0.07 ,0.005+adj_nr,'IPO-','FontWeight','bold','FontSize',20 );
            axestext_c(0.60 ,0.005+adj_nr,'IPO+','FontWeight','bold','FontSize',20 );  
            
         end

end
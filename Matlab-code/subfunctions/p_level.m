%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  p_level.m
%  Daniel Emanuelsson
%  check significance level
%  use with [r,p]=corrcoef_df(detrend(x),detrend(y));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p_string]=p_level(p)

if p>=0.1
    p_string='not_significant';
    
elseif  p<0.1 && p>=0.05
    
    p_string='p<0.1';
    
elseif  p<0.05 && p>=0.01
    
     p_string='p<0.05';
     
elseif  p<0.01 && p>=0.001    
     
     p_string='p<0.01'; 
     
elseif  p<0.001   
         
   p_string='p<0.001';% highly significant 
    
end
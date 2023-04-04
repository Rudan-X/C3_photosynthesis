%
 
%  Copyright   Xin-Guang Zhu, Yu Wang, Donald R. ORT and Stephen P. LONG
%  CAS-MPG Partner Institute for Computational Biology, SIBS/CAS, China 
%  Institute of Genomic Biology and Department of Plant Biology,
%  University of Illinois at Urbana Champaign, United States
%  This file is part of e-photosynthesis. All Rights Reserved.
 
%   e-photosynthesis is distributed for academic research only. 
%   For commercial purpose, please apply for different liscence. 
%   By using this software, you are automatically accepting the academic
%   free licence http://opensource.org/licenses/afl-3.0.php
%   See the content of the licence in the word document: 
%   Academic Free Licence.doc
 
%


%Function [light] = condition; This function is used to store all the required 
%environmental variables, such as light, CO2, O2, humidity as such. This function 
%contains two parts. Part a includes the generic (default) conditions and the 
%second part contains the detailed conditions for different time period.

function fini = Condition (t)

global RUBISCOMETHOD;         %The method for calculation of Rubisco catalyzed reaction
RUBISCOMETHOD = 2;          %1: Use enzyme concentration for calculation
                            %2: Use the michaelis menton and enzyme together for calculation
global VolRatioStCyto
VolRatioStCyto =1; 

                            
%First get the generic conditions
                            
global CO2_cond;
global O2_cond;
global GLight;
global V16;
global Temp_cond;

global Cond_V16;        %This variable is transfered from PSInitial for modificatin of V16, the rate of ATP synthesis. 
CO2Temp = 280;          %CO2 concentation  %ppm
O2Temp = 0.21;          %O2 concentration  %default is 0.21, i.e. 21%. 

CO2_cond = CO2Temp /(3 * 10^4);
O2_cond = O2Temp*1.26;
Temp_cond = 25;

light = 2000;  %light umol m-2 s-1      

%Here the time dependent variable is regulated. 
global tglobal;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  Here define how many interval needed for the experiments  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NumberInterval = 10;

global NumInter_draw; 
NumInter_draw = NumberInterval;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This is an experimental protocol for doing repeatative experiment
%
%Tinter = tglobal/NumberInterval;
%
%FirstMet = 0;
%
%    for index = 1:NumberInterval 
%            b = index * Tinter;
%        if t <= b & FirstMet == 0
%            
%            modifier = 1 * index; 
%            
%            global BF2XanCycle_pHl;
%
%            global Xan2Stom_ABA;
%
%                
%            %Light regulation
%%            light =  100 * index; 
%            light = 2000; 
%            
%            %CO2 regulation
%%            temp = 280;
%            temp = 1000-100 * (index-1); 
%            
%            %O2 regualtion
%            O2Temp = 0.21 ;   
%            
%            CO2_cond = temp /(3 * 10 ^ 4);   
%            O2_cond = O2Temp*1.26;
%            
%            %Regulation of Vmax
%            %Global RuACT_RC;
%            %RuACT_RC(9) = 25 + 25  * index/2;
%            
%            FirstMet = 1;                
%            
%        end
%    end
%        
%        
%
%%%%PAM MEASUREMENT
%%if t < 20
%%    light = 7;
%%else
%%    light=500;
%%end
%%light = 500;
%%%%PAM MEASUREMENT
%%  StepL = 20; 
%%    if t<1
%%        light = 8000;
%%    elseif t>1 & t < 20
%%        light = 1;
%%    else 
%%        index = floor(t/StepL);
%%        if t > (StepL*index) & t < (StepL*index+1)
%%            light = 8000; 
%%        else 
%%            light = 500; 
%%        end
%%    end   
%
%%
%%%Experiment 4 in Zhu 2012

%if t<200
%    light = 1000;
%elseif t>200 & t<400
%    light = 100;
%else
%    light = 1000;
%end

% if t<200
%     CO2_cond = 1000;
% elseif t>200 & t<400
%     CO2_cond = 900;
% elseif t>400 & t<600
%     CO2_cond = 800;
% elseif t>600 & t<800
%     CO2_cond = 700;
% elseif t>800 & t<1000
%     CO2_cond = 600;
% elseif t>1000 & t<1200
%     CO2_cond = 500;
% elseif t>1200 & t<1400
%     CO2_cond = 400;
% elseif t>1400 & t<1600
%     CO2_cond = 300;
% elseif t>1600 & t<1800
%     CO2_cond = 200;
% else
%     CO2_cond = 100;
% end

t_inter_ca=200: 200: 1800;

t_Ca=[50:50:200, 350, 500, 750, 1000, 1500];

ind=find(t_inter_ca>=t,1,'first');
% fprintf("Time: %4.2f, Ca: %4.2f\n",t, t_Ca(ind))
CO2_cond=t_Ca(ind);

CO2_cond = CO2_cond/(3 * 10^4);
%%%
%disp([t CO2_cond])
GLight   = light;
fini = 1;   
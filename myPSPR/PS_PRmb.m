%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
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
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function PS_PR_DYDT = PS_PRmb(t,PS_PRs,PS_PR_Param)

PSs = zeros(15,1);
PrS = zeros(13,1);

for m=1:4
    PSs(m)=PS_PRs(m);
end

for m = 5:14
    PSs(m+1)=PS_PRs(m);
end

PSs(5) = PS_PRs(24);

for m = 15:16
    PrS(m-14) = PS_PRs(m);
end
PrS(3) = PS_PRs(2);            

for m = 17:23
    PrS(m-13)=PS_PRs(m);
end


PrS(11) = PS_PRs(1);            %RUBP
PrS(12) = PS_PRs(11);           %CO2
PrS(13) = PS_PRs(12);           %O2
    
PR2PS_Pgca = PrS(4);            %FOr transfering information between PR to PS.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%2. Add exprimental conditions here; Conditions like light, temperature, CO2, O2 concentration should be added here %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fini = Condition (t);


PR_v = zeros(10,1);
PS_v = zeros(17,1);

PS_Param = zeros(2,1);
PS_Param (1) = PS_PR_Param;
PS_Param(2) = PR2PS_Pgca;

PS_v = PSRate(t, PSs,PS_Param);

PR_Param = zeros(2,1);
PR_Param(1) = PS_PR_Param;                  %To indicate that the calcualtion is using the combined model
                                 %for the PS-PR combined model. 0: Combined model; 1: Separate model
global PS2PR_Pi;
PR_Param(2) = PS2PR_Pi;

PR_v = PRrate(t,PrS,PR_Param);


%Assign the rate of reaction that is calculated from the photosynthesis and photorespiration routine.

v1	=	PS_v(1)	;
v2	=	PS_v(2)	;
v3	=	PS_v(3)	;
NONE=	PS_v(4)	;
v5	=	PS_v(5)	;
v6	=	PS_v(6)	;
v7	=	PS_v(7)	;
v8	=	PS_v(8)	;
v9	=	PS_v(9)	;
v10	=	PS_v(10);
v13	=	PS_v(11);
v16	=	PS_v(12);
v23	=	PS_v(13);
v31	=	PS_v(14);
v32	=	PS_v(15);
v33	=	PS_v(16);
v24 =   PS_v(17);
v25 =   PS_v (18); 

v111 = PR_v(1)  ;
v112 = PR_v(2)  ;
v113 = PR_v(3)  ;
v121 = PR_v(4)  ;
v122 = PR_v(5)  ;
v123 = PR_v(6)  ;
v124 = PR_v(7)  ;
v131 = PR_v(8)  ;
v1in = PR_v(9)  ;
v2out= PR_v(10) ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%5. This part exchange informations from two systems.%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global PSPR2RA_v1;
global PSPR2RA_v13;
global PSPR2RA_v111;

PSPR2RA_v1 = v1;
PSPR2RA_v13 = v13;
PSPR2RA_v111 = v111;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%6.  Calculation of the mass balance equations %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$

tmp = zeros(23,1);


tmp(1) = v13-v1-v111;


tmp(2) = 2*v1 - v2 - v32 + v113 + v111;


tmp(3) = v2 - v3;


tmp(4) = v3 - 2* v5 - v7 - v8 - v10 - v31 - v33; 


tmp(5) = v5-v6;    


tmp(6) = v7-v8;


tmp(7) = v9-v10;


tmp(8) = v8 - v9;


tmp(9) = v16 - v2 - v23 - v13- v113 - v25;    
  


tmp(10) = 0;


tmp(11) = 0;


tmp(12) = 0;


tmp(13) = v6 - v7 - v23 + v25;


tmp(14) = v7 + v10 * 2 - v13;


tmp(15) = v1in  - v113;


tmp(16) =  v112 - v2out;


tmp(17) = v111 - v112;


tmp(18) = v2out - v121;


tmp(19) = v121 - v122- v124;


tmp(20) = v131 - v122;


tmp(21) = v122 + v124 - 2*v131;

tmp(22) = v122 - v123;
tmp(23) = v123 - v1in;

tmp(24) = v23 - v24;   

PS_PR_DYDT = tmp;
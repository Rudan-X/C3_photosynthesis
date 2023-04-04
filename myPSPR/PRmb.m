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
%   See the content of the licence  in the word document: 
%   Academic Free Licence.doc
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




function PR_DYDT = PRmb(t,PrS,PR_Param)

fini = Condition (t);

Velocity = zeros(10,1);
Velocity = PRrate(t,PrS,PR_Param);

v111 = Velocity(1)  ;
v112 = Velocity(2)  ;
v113 = Velocity(3)  ;
v121 = Velocity(4)  ;
v122 = Velocity(5)  ;
v123 = Velocity(6)  ;
v124 = Velocity(7)  ;
v131 = Velocity(8)  ;
v1in = Velocity(9)  ;
v2out = Velocity(10)  ;


tmp = zeros(11,1);


tmp(1) = v1in  - v113;


tmp(2) =  v112 - v2out ;


tmp(3) = v111+ v113 - 0.5;

tmp(3) = 0;


tmp(4) = v111-v112;


tmp(5) = v2out - v121;



tmp(6) = v121 - v122- v124;


tmp(7) = v131 - v122;


tmp(8) = v122 + v124 - 2*v131;


tmp(9) = v122 - v123;


tmp(10) = v123 - v1in;


tmp(11) = 0.3 - v111;


tmp(12) = 0;


tmp(13) = 0;


PR_DYDT = tmp;

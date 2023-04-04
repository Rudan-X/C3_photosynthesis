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

function ModelComb = IniModelCom;

global RuACT_EPS_com;
RuACT_EPS_com = 0;

global BF_FI_com;            %The combination of BF and FI model 
BF_FI_com = 0;

global PR_PS_com;    %This is a variable indicating whether the PR model is actually need to be combined with PS or not. If 1 then means combined; 0 means not. 
PR_PS_com = 0;

global FIBF_PSPR_com; %1 means that the overall EPS model is used. 0 means partial model of FIBF is used. 
FIBF_PSPR_com = 0;    

global ATPActive;
ATPActive = 0;

global RedoxReg_RA_com;
RedoxReg_RA_com = 0;

global XanCycle_BF_com;
XanCycle_BF_com = 0;

global RROEA_EPS_com;
RROEA_EPS_com = 0;

global StomCond_TrDynaPS_com;
StomCond_TrDynaPS_com = 0;

global PSPR_SUCS_com;
PSPR_SUCS_com = 0;  

global trDynaPS_SUCS_com;
trDynaPS_SUCS_com = 0;

global EPS_SUCS_com;
EPS_SUCS_com = 0;

ModelComb = 1;
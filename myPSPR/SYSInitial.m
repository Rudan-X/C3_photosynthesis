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

%SYSInitial.m
%This routine initialze the parameters used for all the routines.

function fin = SYSInitial(Begin)

%The total runing time is set as tglobal
global tglobal;
global Primer; 
Primer = 300; 

tglobal = 600; 

 
global options1
options1 = odeset('RelTol',1e-4);

global PS12ratio;        
PS12ratio = 1;

global AVR;         
AVR = 30; 

global GP; 
GP = 0; 

fin = 1; 
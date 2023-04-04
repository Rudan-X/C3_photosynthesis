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



function PSs = PSInitial(my_params)

global PS_C_CA;             %  Global constant for the total adenylates
global PS_C_CP;             %  Global constant for the total phosphate
global PS_C_CN;             %  Global constant for the total NADP+NADPH
global PS_PEXT;             %  Global constant for the cytosolic Phosphate concentration;

PS_C_CP	= 15; 
PS_C_CA	=1.5;        
PS_C_CN	=1;
PS_PEXT	=0.5;

global PSPR_RA_CA;
PSPR_RA_CA = PS_C_CA;

RuBP	= 2.000;
PGA	    = 2.400;
DPGA	=0.0011;
T3P	    = 0.5;
ADPG	= 0.005;          
FBP	    =0.670;
E4P	    =0.050;
S7P	    = 2.000;        
SBP	    = 0.300;
ATP	    = 0.68;             
NADPH	= 0.21;         
HexP    = 2.2;          
PenP    = 0.25;

CO2	    = 0.012;          
O2	    = 0.21*1.26 ;        

PSs(1)	=RuBP;
PSs(2)	=PGA;
PSs(3)	=DPGA;
PSs(4)	=T3P;
PSs(5)	=ADPG;
PSs(6)	=FBP;
PSs(7)	=E4P;
PSs(8)	=S7P;
PSs(9)	=SBP;
PSs(10)	=ATP;
PSs(11)	=NADPH;
PSs(12)	= CO2;
PSs(13)	= O2;
PSs(14) = HexP;
PSs(15) = PenP;


%Initialize the constants for the different reactions

global	KM11	;
global	KM12	;
global	KM13	;
global  KI11    ;
global  KI12    ;
global  KI13    ;
global  KI14    ;
global  KI15    ;

global	KM21	;
global	KM22	;
global  KM23    ;

global	KM31a	;
global	KM32b	;

global	KM41	;
global	KM42	;
global  KE4     ;

global	KM51	;
global	KM52	;
global	KM53	;
global  KE5     ;

global	KM61	;
global  KI61    ;
global  KI62    ;

global	KM71	;
global  KM72    ;
global	KM73	;
global  KM74    ;

global	KM8	    ;
global  KM81    ;
global  KM82    ;

global	KM9	    ;
global  KI9     ;

global	KM10	;
global  KM101   ;
global  KM102   ;
global  KM103   ;

global	KE11	;

global	KE12	;

global	KM131	;
global	KM132	;
global	KI131	;
global	KI132	;
global	KI133	;
global	KI134	;
global	KI135	;

global	KM161	;
global	KM162	;

global	KE21	;

global	KE22	;



global	KM311	;
global	KM312	;
global	KM313	;

global	KM32	;
global	KM33	;


global KE6;
global KE7;
global KE8;
global KE9;
global KE10;
global KE13;
global KE16;

global KM103;
global KM163;


%% my code
% KM11	=	0.0115;		%	CO2	1	RuBP+CO2->2PGA 
% KM12	=	0.222;		%	O2	1	RuBP+CO2->2PGA 
% KM13	=	0.02;		%	RuBP	1	RuBP+CO2->2PGA 
KM11=my_params(1);
KM12=my_params(2);
KM13=my_params(3);
%%
KI11    =   0.84   ;    %PGA
KI12    =   0.04   ;       %FBP
KI13    = 0.075 ;       %SBP
KI14    = 0.9   ;       %Pi
KI15    = 0.07  ;       %NADPH

KM21	=	0.240;		%	PGA	2	PGA+ATP <-> ADP + DPGA
KM22	=	0.390;		%	ATP	2	PGA+ATP <-> ADP + DPGA
KM23    =   0.23  ;        % ADP     

KM31a	=	0.004;		%	BPGA	3	DPGA+NADPH <->GAP + OP+NADP 
KM32b	=	0.1	;	    %	NADPH	3	DPGA+NADPH <->GAP + OP+NADP

KM41	=	2.5	;	    %	DHAP	4	DHAP <->GAP 
KM42	=	0.68;		%	GAP	4	DHAP <->GAP
KE4     =   1/0.05;       %  Using the value from Patterson

KM51	=	0.3	;	    %	GAP	5	GAP+DHAP <->FBP
KM52	=	0.4	;	    %	DHAP	5	GAP+DHAP <->FBP
KM53	=	0.02;		%	FBP	5	GAP+DHAP <->FBP     %Original Value: 0.02
KE5     =   7.100;          %Defult: 7.1

KM61	=	0.033;		%	FBP	6	FBP<->F6P+OP
KI61    = 0.7   ;       %  F6P       
KI62    = 12    ;       %  Pi 
KE6     =   6.66 * 10^5;    %The equilibrium constant for this reaction        %New    mM     Laisk or Bassham and Krause 1969 BBA

KM71	=	0.100;		%	Xu5P	7	F6P+GAP<->E4P+Xu5P      %jn
KM72	=	0.100;		%	E4P	7	F6P+GAP<->E4P+Xu5P
KM73    = 0.1;         %  F6P This value was based on estimate
KM74    = 0.1000;         %Estimate for GAP ORIGINAL 0.1 
KE7     =   0.1 ;       %The equilibrium constant for this reaction             %New           Laisk  Bassham and Krause 1969 BBA

KM8	    =	0.02;		%	SBP	8	E4P+DHAP<->SBP
KM81    = 0.4   ;       %DHAP
KM82    = 0.2   ;       %E4P estimate
KE8     = 1.017 ;     %The equilibrium constant for this reaction                  %New    mM-1         Laisk  Bassham and Krause 1969 BBA. Default: 1.107

KM9	    =	0.05;		%	SBP	9	SBP<->S7P+OP    
KI9     = 12    ;       %  The inibintion constant for Pi; 
KE9     =   6.66 * 10^5 ; %The equilibrium constant of this reaction           %New   mM      Laisk  Bassham and Krause 1969 BBA

KM10	=	0.5	;	    %	R5P	10	S7P+GAP<->Ri5P+Xu5P
KM101   =   0.1 ;       %  Xu5P
KM102   = 0.09 ;       %  Estimate for GAP
KM103   = 0.015 ;        %  Estimate for S7P                                    %New 
KE10    = 1/0.85 ;      %  The equilibrium constant for this reaction          %New From Laisk or Bassham and Krause 1969 BBA

KE11	=	0.4	;	    %	Equilibrium Constant	11	Ri5P<-->Ru5P
KE12	=	0.67;		%	Equilibrium Constant	12	Xu5P<-->Ru5P

% KM131	=	0.05;		    %	Ru5P	13	Ru5P+ATP<->RuBP+ADP
% KM132	=	0.059;		    %	ATP	13	Ru5P+ATP<->RuBP+ADP
% KI131	=	2	;			%	PGA	13	Ru5P+ATP<->RuBP+ADP
% KI132	=	0.7	;			%	RuBP	13	Ru5P+ATP<->RuBP+ADP
% KI133	=	4	;			%	Pi	13	Ru5P+ATP<->RuBP+ADP
% KI134	=	2.5	;			%	ADP	13	Ru5P+ATP<->RuBP+ADP

%% my code
KM131	=	my_params(4);		    %	Ru5P	13	Ru5P+ATP<->RuBP+ADP
KM132	=	my_params(5);		    %	ATP	13	Ru5P+ATP<->RuBP+ADP
KI131	=	my_params(6);			%	PGA	13	Ru5P+ATP<->RuBP+ADP
KI132	=	my_params(7);			%	RuBP	13	Ru5P+ATP<->RuBP+ADP
KI133	=	my_params(8);			%	Pi	13	Ru5P+ATP<->RuBP+ADP
KI134	=	my_params(9);			%	ADP	13	Ru5P+ATP<->RuBP+ADP

%%
KI135	=	0.4	;			%	ADP	13	Ru5P+ATP<->RuBP+ADP
KE13    =   6.846 * 10^3;   %  The equilibrium constant for this reaction  %New From Laisk or Bassham and Krause 1969 BBA

KM161	=	0.014;		%	ADP	16	ADP+Pi<->ATP
KM162	=	0.3;		%	Pi	16	ADP+Pi<-> ATP
KM163   =   0.3;        %  ATP 16  ADP+Pi<-> ATP                           %New       Based on Laisk  
KE16    =   5.734;      %  The equilibrium constant for this reaction      %NEW, From Laisk or Bassham and Krause 1969 BBA


KE21	=	2.3;		%	Equilibrium constant	21	F6P<->G6P 
KE22	=	0.058;		%	Equilibrium constant	22	G6P<->G1P

%KM231	=	0.08;		%	G1P	23	G1P+ATP+Gn<->PPi+ADP+Gn+1
%KM232	=	0.08;		%	ATP	23	G1P+ATP+Gn<->PPi+ADP+Gn+1
%KA231	=	0.1;		%	PGA	23	G1P+ATP+Gn<->PPi+ADP+Gn+1
%KA232	=	0.02;		%	F6P	23	G1P+ATP+Gn<->PPi+ADP+Gn+1
%KA233	=	0.02;		%	FBP	23	G1P+ATP+Gn<->PPi+ADP+Gn+1
%KI23	=	10;		    %	ADP	23	G1P+ATP+Gn<->PPi+ADP+Gn+1

KM311	=	0.077;		%	DHAP	31	DHAPi<->DHAPo
KM312	=	0.63;		%	Pi	31	DHAPi<->DHAPo
KM313	=	0.74;		%	Pext	31	DHAPi<->DHAPo
KM32	=	0.25;		%	PGA	32	PGAi<->PGAo
KM33	=	0.075;		%	GAP	33	GAPi<->GAPo


%Now put in the constant for the new ADPG Pyrophosphorylase and starch
%synthase

%ATP + Glucose-1-Phosphate --> ADPG + PPi

global	KM231	;
global	KM232	;
global	KM233	;
global	KM234	;
global	KA231	;
global	KI231	;
global  KVmo    ;
global  KE23    ;

KM231	=	0.031;		%	G1P	23	G1P+ATP<->ADPG + PPi        Laisk et al 1989
KM232	=	0.045;		%	ATP	23	G1P+ATP<->ADPG + PPi        Laisk et al 1989
KM233	=	0.14;		%	ADPG	23	G1P+ATP<->ADPG + PPi        Laisk et al 1989
KM234	=	0.8;		%	PPi	23	G1P+ATP<->ADPG + PPi        Laisk et al 1989
KE23	=	7.6 * 10^(-3);		%	PPi	23	G1P+ATP<->ADPG + PPi        Laisk et al 1989


KA231	=	0.23;		%	PGA	23	G1P+ATP<->ADPG + PPi        Laisk et al 1989
KI231   =   0.9 ;       %  Pi	23	G1P+ATP<->ADPG + PPi        Laisk et al 1989
KVmo    =   0.007;      %  The minimum maximum velocity        Laisk et al 1989
%ADPG --> ADP + Gn     %The starch synthesis reaction 24.     Laisk et al
%1989

global KM241;
global KM242;
global KE24;

KM241   =   0.2 ;   %  ADPG    ADPG --> ADP + Gn       Laisk et al 1989
KM242   =   0.6 ;   %  ADP     ADPG --> ADP + Gn       Laisk et al 1989
KE24    =   7.4 * 10^5  ;   %  ADP     ADPG --> ADP + Gn       Laisk et al 1989

global KE25; 
KE25 = 1.2 * 107;

%Initialize the Vmax for different reactions

global	V1	;
global	V2	;
global	V3	;
global	V5	;
global	V6	;
global	V7	;
global	V8	;
global	V9	;
global	V10	;
global	V11	;
global	V12	;
global	V13	;
global	V16	;
global	V21	;
global	V22	;
global	V23	;
global	V31	;
global	V32	;
global	V33	;
global  V24  ;


global GP; 
if GP ==0
	%FC is a fussl factor here. 
	FC = 1;       %Defulat is 2.5. 
	fc16 = 1;     %default 1.5.
	%Initialize the values of the global variables
	SC = 1;        %Scalling coefficient for the stroma volume per mg chl. defualt 2
	SC1 = 1; 
	
    SC222 = 2; 
	STOM1  = 1; 
	STOM2 = 1;
	
	V1		=	 2.93 * SC1 /STOM1 ;	%	(Harris & Koniger, 1997)	
	V2		=	 30.15 * SC	* STOM2 ;	%	(Harris & Koniger, 1997)	
	V3		=	 4.04 * SC * STOM2;%1.57*SC	    ;	%	(Harris & Koniger, 1997)
	V5		=	1.22*SC	    ;	%	(Harris & Koniger, 1997)	
	V6		=	0.734*SC/STOM1	;	%	(Harris & Koniger, 1997)	
	V7		=	3.12*SC * 4;	%	(Harris & Koniger, 1997)	
	V8		=	1.22*SC	   ;	%	(Harris & Koniger, 1997)	
    V9		=	0.32 *3*0.5  ;       %0.17*SC *FC	;	%	(Harris & Koniger, 1997) *3. 
	V10		=	V7;	%	(Harris & Koniger, 1997)	
    V13		=	10.81*SC1	;	%	(Harris & Koniger, 1997)
	V16		=	5.47;       %(Aflalo & Shavit, 1983, Davenport & McLeod, 1986)	
	V23		=	2  ;
    
    
%	V1	=	2.913930914	;
%	V2	=	30.1408176	;
%	V3	=	4.039482839	;
%	V5	=	1.218890457	;
%	V6	=	0.726259575	;
%	V7	=	3.122215944	;
%	V8	=	1.218890364	;
%	V9	=	0.324190582	;
%	V10	=	V7;
%	V13	=	10.83475664	;
%	V23	=	0.26684349	;
%  V16 = 5; 

end
    V24     =   2;           
	V31		=	1.0  ;      
	V32		=	1.0 ;       
	V33		=	1.0 ;        

    
global Cond_V16;        %This parameter is used to modifying the V16 if needed. 
Cond_V16 = V16; 

global PS2SUCSV32;
PS2SUCSV32 = V32;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Here is the location for transfering variables%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global PS2RA_RuBP_ini;
PS2RA_RuBP_ini = RuBP;

global BF_FI_com;
global PS2BF_ATP;
global PS2BF_ADP;
global PS2BF_Pi;

PS2BF_ATP = ATP;
PS2BF_ADP = PS_C_CA - ATP;

global PS2PR_V1;
PS2PR_V1 = V1;

global CO2_PS2StomCond; 
CO2_PS2StomCond = CO2; 

global PS2SUCS_PGA;
PS2SUCS_PGA = PGA; 
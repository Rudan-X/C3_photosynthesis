function z = optim_obj(Variable,model,error_type)
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


Begin = 1;
fin = SYSInitial(Begin);
global options1;
global tglobal;
time = 1800;% tglobal;


global PS_PR_OLDTIME;
global PS_PR_TIME_N;
global PS_PR_VEL;

PS_PR_OLDTIME = 0;
PS_PR_TIME_N = 1;
PS_PR_VEL = zeros(27,1);        

global PS_OLD_TIME;
global PS_TIME_N;
global PS_VEL;
PS_OLD_TIME = 0;
PS_TIME_N= 0;
PS_VEL = zeros(1,1);

global PR_OLD_TIME;
global PR_TIME_N;
global PR_VEL;
PR_OLD_TIME = 0;
PR_TIME_N = 1;
PR_VEL = zeros(1,1);

%%%%%%%%%%%%%%%%%%%%%%%%
%  Initialation step %
%%%%%%%%%%%%%%%%%%%%%%%%
Begin = 1;

ModelComb = IniModelCom;       

global PR_PS_com;     
PR_PS_com = 1;

PS_PRs = PS_PRIni(Variable);

PS_PR_Param = 0;

global CO2A;
CO2A = zeros(5,1);
tspan=[0,1800];
ode_sol = ode15s(@PS_PRmb,tspan,PS_PRs,options1,PS_PR_Param);

if ~isempty(ode_sol) && max(ode_sol.x)==max(tspan)
    global AVR; 
    Coeff = AVR ;
    
    CarbonRate = PS_VEL(2,:)*Coeff;
    CO2Release = PR_VEL(:,9)*Coeff;
    
    NetCarbon = CarbonRate' - CO2Release;
    
    % plot(PS_VEL(1,:),NetCarbon,'.');
    time=PS_VEL(1,:);
    V=[200: 200: 1800]'-5;
    N=time';
    A = repmat(N,[1 length(V)]);
    [~,closestIndex] = min(abs(A-V'));
    simulated=NetCarbon(closestIndex);

    measured=model.ACIC(:,2);
    sd=model.ACIC(:,3);
    
    if strcmp(error_type,"relative_error")
        z=sum(abs(measured-simulated)./measured,'omitnan');
    elseif strcmp(error_type,"chi_square")
        z=sum((measured-simulated).^2./sd.^2,'omitnan');
    end


else
    z=1e10;
end

% calculate carbon assimilation


fprintf('%s: %d\n',error_type,z)


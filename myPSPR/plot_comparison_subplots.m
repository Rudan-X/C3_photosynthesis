function z = plot_comparison_subplots(Variable,acc_i,title_name)

%%
T = readtable("aci.xlsx");
ind_table=find(T{:,"Genotype"}==acc_i);

fprintf("Optimizing for accession: %s\n",acc_i)

ACIc=table2array(T(ind_table,[2:3,8]));


model.ACIC=ACIc;

%% Objective function
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
    
    z=sum((measured-simulated).^2./sd.^2,'omitnan');

else
    z=1e10;
end


%%%%%%%%%%%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%%%%%

x=[50:50:200, 350, 500, 750, 1000, 1500]; % Ca
yu=measured+sd;
yl=measured-sd;
x=x';% x must be a column
yu=yu';
yl=yl'; % yl AND yu must be rows.
xfill=[x',flipud(x)']; % xfill must be a single row: the second curve is reversed. It should be like [1 2 3 3 2 1] It fills between two curves. It must be in the same row, and
yfill=[yl, fliplr(yu)];


fill(xfill,yfill,[0 0.4470 0.7410],'LineStyle','none','FaceAlpha',.5);
hold on;
plot(x,measured,'-o','color',[0 0.4470 0.7410],'MarkerSize',8, 'linewidth',1.5)
hold on 
plot(x,simulated,'*','color',[0.9290 0.6940 0.1250],'MarkerSize',8, 'linewidth',1.5)
% ylabel("A (µmol/m2/s)")
% xlabel("time (min)")
% acc_lab=replace(acc_i,"_","-");
% title(strcat(acc_lab,", error: ",string(z))) 
title(title_name,'FontSize', 17)
ylim([-1.6, 26]);

fprintf('Chi square: %d\n',z)



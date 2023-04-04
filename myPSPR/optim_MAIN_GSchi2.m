function []= optim_MAIN_GSchi2(arg_ind)

%% Load ACI curve
T = readtable("aci.xlsx");
accessions=string(unique(T{:,"Genotype"}));
% accessions=["c24","col"];

global acc_i
acc_i=accessions(arg_ind);
ind_table=find(T{:,"Genotype"}==acc_i);

fprintf("Optimizing for accession: %s\n",acc_i)

ACIc=table2array(T(ind_table,[2:3,8]));


model.ACIC=ACIc;

init_sol=[0.0115,0.222,0.02, 0.05, 0.059, 2, 0.7, 4, 2.5, 0.4]	;


fun = @(x)optim_obj(x,model,"chi_square");

%% Set up shared variables with outfun

display("Global search")
lb=1e-3*ones(length(init_sol),1);
ub=1e2*ones(length(init_sol),1);

gs = GlobalSearch('MaxTime',60*60*24*7); % 21 days
options = optimoptions('fmincon','algorithm','interior-point','OutputFcn',@optim_output_GS,'Display','iter');

problem = createOptimProblem('fmincon','objective',fun, 'x0', init_sol,...
    'lb',lb, 'ub',ub,'options',options); 



tic
[xsol,fval,exitflag,output] = run(gs,problem);
toc

filen=strcat("optim_GSresult_",acc_i,".mat");
save(filen,"xsol","fval","exitflag","output")

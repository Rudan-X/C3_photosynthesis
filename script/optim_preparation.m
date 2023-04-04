clear
clc

%% Model construction
model_building()

% parse_rules.m adds new fields to model: .grRules, .genes(from araCORE),
% .protein_abun extracted from Piques's papaer (doi:10.1038/msb.2009.68),
% .enzymes which catalyze each reaction.
model=model_parse_rules(model); 

% decompose_model.m returns the decomposed model with their corresponding
% fields, but .enzymes and .protein_abun are kept in their original scale
[model_d,model0]=model_decomposing(model);

model_d=convertToIrreversible(model_d);

%% Getting experimental data

% 1. Load experimental data: in units mili mol/L
% number of mol per the amount unit requested:
mol_per_amount_unit=1e-3; % 0.001 mili mol per mol: returned amount units= mili mol

% returns: model_d.exp_meta=[concentration, std dev]
[model_d]=load_experiment_data_total(model_d,mol_per_amount_unit); 

% 2. Load light irradiance in units: mili mol/L/h
% number of seconds per the time unit requested:
second_per_time_unit=3600; % 3600 second per hour: returned time units=hour

% load_irradiance returns model_d.constant, a structure storing the fixed flux for the reaction
% which import the photons
model_d=load_irradiance(model_d,mol_per_amount_unit,second_per_time_unit);

% 3.Load constant variables during ODE
% Record the metabolites with fixed concentrations, e.g. CO2, Pi, O2, which
% are not simulated in ODE, but set to constant values, in units mili mol/L
% returns model_d.constant=[day,night]

% [model_d]=load_ODE_cte(model_d,mol_per_amount_unit); 


% find the index of irradiance at time point 0 for condition 1
model_d.cond=1;
% start time: before irradiation
light_influx0=model_d.irradiance{model_d.cond}(1,1);


% 4. Collect initial concentrations for some metabolites
% collect_concentration.m returns model_d.metslb and model_d.metsub which
% record the concentration for some metabolites

% model_d=load_collected_concentration(model_d);


% 5. Building the fields needed for optimization

start_time=-0.5; % starting time point
model_d.tf=24; % final time point

model_d.cond=1; % numbers of conditions evaluated
model_d.k0=length(model_d.mets); % position where k starts


model_d.c3d_meas=cell(length(model_d.cond),1);
model_d.sd3d=cell(length(model_d.cond),1);

i=0;
model_d.irra_ode=cell(length(model_d.cond),1);
for cond=1:length(model_d.cond)
    if sum(cond==model_d.cond)~=0
        i=i+1;
        index=1:size(model_d.exp_meta_total{cond,1},1);
        
        model_d.c3d_meas{i}=model_d.exp_meta_total{cond,1}(index,1:(end-1))';
        
        model_d.sd3d{i}=model_d.exp_meta_total{cond,2}(index,1:(end-1))';
        model_d.irra_ode{i}=model_d.irradiance{cond}(index,:);
    end
end

lb_k=1e-4*ones(length(model_d.rxns),1);
ub_k=1e3*ones(length(model_d.rxns),1);


% lower and upper bounds for the variable
lb=[1e-4*ones(length(model_d.mets),1);lb_k]; 
ub=[5000*ones(length(model_d.mets),1);ub_k];

% 6. Optimization constraints
[model_d,Aeq,beq,count_met]=load_constraint(model_d);


A=[];
b=[];

% 7. Get some initial values for our variable
nvars=length(lb);
prob=struct();
prob.lb=lb;
prob.ub=ub;
prob.obj=zeros(nvars,1);

prob.A=[sparse(Aeq); sparse(A)];
prob.rhs = full([beq;b]);
prob.vtype = repmat('C', nvars, 1);
prob.sense = [repmat('=',size(Aeq,1),1);repmat('<',size(A,1),1)];

sol=gurobi(prob);

%% 8. Change some of initial values to the ones that showed starch accumulation

default_m=0.5;
default_more=5;
default_k=1000;


[x0,ks,model_d]=change_init_sol(model_d,sol,default_m,default_more,default_k,count_met); 
% change_init_sol.m adds model_d.irra_value as well

model_d.fix_ind=1;
model_d.sim_ind=setdiff(1:length(model_d.mets),1);


% changed lb and ub for Ks and initial light irradiance
lb=[1e-4*ones(length(model_d.mets),1);1e-2*ones(length(model_d.rxns),1)]; 
ub=[1000*ones(length(model_d.mets),1);3e4*ones(length(model_d.rxns),1)];



save("./server_script/prepared_data.mat")
%% check how initial variables work

% model_d.fix_ind=1;
% model_d.sim_ind=setdiff(1:length(model_d.mets),1);
% 
% x0=x0(model_d.sim_ind);
% 
% 
% tspan=model_d.irra_ode{model_d.cond}(:,2);
% model_d.irradiance=model_d.irradiance;
% 
% startT=datetime('now');
% overfunc=@(sim_t,sim_x)ode_stopfunc(sim_t,sim_x,startT);
% options=odeset('NonNegative',1:length(model_d.sim_ind),'RelTol',1e-5,'AbsTol',1e-5,'Events',overfunc); %  
% tic
% ode_sol=ode15s(@(sim_t,sim_x)ode_rates(sim_t,sim_x,ks,model_d),tspan,x0,options); 
% toc
% 
% 
% maxtime=max(ode_sol.x);
% ind=find(tspan<=maxtime,1,'last');
% pred_concent=zeros(length(model_d.mets),length(1:ind));
% pred_concent(model_d.sim_ind,:)=deval(ode_sol,tspan(1:ind));
% 
% irra=model_d.irra_ode{1}(:,1);
% pred_concent(1,irra~=0)=model_d.irra_value;
% 
% plot_comparison(pred_concent,model_d,"all")



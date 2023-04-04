T = readtable("aci.xlsx");
accessions=string(unique(T{:,"Genotype"}));
% accessions=["c24","col"];

%%
fig=figure;
p=0;
title_names=["C24","Col-0"];
for k=2:3
    p=p+1;
    my_sub=subplot(1,2,p);
    acc_i=accessions(k);
    
    filen=strcat("optim_GSresult_",acc_i,".mat");
    res=load(filen);
    optimized=res.xsol;
    plot_comparison_subplots(optimized,acc_i,title_names(p)) 
    set(my_sub, 'fontsize', 12); %,'FontWeight','bold'
    fprintf("Objective value: %4.3f\n", res.fval)
end

han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
% ylabel(han,"",'Interpreter','Latex','FontSize', 17,'String','{\boldmath$A [\mu mol \, m^{-2}\, s^{-1}]$}');
% xlabel(han,"",'Interpreter','Latex','FontSize', 17,'String','{\boldmath$C_i [ppm]$}');
ylabel(han,"A (µmol/m2/s)",'FontSize', 15);
xlabel(han,"Ca (µmol/mol)",'FontSize', 15);

filen=strcat("./report/C3results_2only.png");
% set(gcf, 'PaperPosition', [0 0 20 7]); %
set(gcf, 'PaperPosition', [0 0 12 5]);

saveas(gcf,filen)




%%
k=1;
figure
acc_i=accessions(k);

filen=strcat("optim_GSresult_",acc_i,".mat");
res=load(filen);
optimized=res.xsol;
plot_comparison_subplots(optimized,acc_i) 
legend(["measured +/- std","measured","simulated"],"FontSize",17,'FontWeight','bold')


filen=strcat("./report/legend.png");
set(gcf, 'PaperPosition', [0 0 20 7]); %
saveas(gcf,filen)

%%

my_var=zeros(9,10);
max_A=zeros(1,10);
for k=1:10
    acc_i=accessions(k);
    filen=strcat("optim_GSresult_",acc_i,".mat");
    res=load(filen);
    my_var(:,k)=res.xsol(1:9);

    ind_table=find(T{:,"Genotype"}==acc_i);
    ACIc=table2array(T(ind_table,[2:3,8])); 
    max_A(k)=max(ACIc(:,2));

end

kine_names=["RuBisCO KM(CO2)","RuBisCO KM(O2)","RuBisCO KM(RuBP)","PRK KM(Ru5P)",...
    "PRK KM(ATP)","PRK KM(PGA)","PRK KM(RuBP)","PRK KM(Pi)","PRK KM(ADP)"];

fig=figure;
ind_new=[2,3];
init_val=[0.0115,0.222,0.02, 0.05, 0.059, 2, 0.7, 4, 2.5];
for var=1:9
    my_sub=subplot(3,3,var);
    Y=my_var(var,ind_new);
    bar(Y)
    text(1:length(Y),Y,num2str(round(Y',3)),'vert','bottom','horiz','center'); 
    set(my_sub, 'fontsize', 12); %,'FontWeight','bold'
    title(strcat(kine_names(var),": ",string(init_val(var))),'FontSize', 12)
    ylim([0,1.3*max(Y)])

    xticklabels(["C24","Col-0"])
%     xtickangle(45)
end

han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
% xlabel(han,"",'Interpreter','Latex','FontSize', 15,'String','{\boldmath$Time \, points$}');
% xlabel(han,"Accessions",'FontSize', 15);
ylabel(han,"KM (mM)",'FontSize', 15);
xlabel(han,"Accessions",'FontSize', 15);

filen=strcat("./report/KMdistribution_2only_initial.png");
set(gcf, 'PaperPosition', [0 0 10 9]); %
saveas(gcf,filen)


%% Correlation between the maximal A and the KMs

figure
for var=1:9
    my_sub=subplot(3,3,var);
    plot(max_A,my_var(var,:),'o')
    % title(strcat(acc_lab,", error: ",string(z))) 
    set(my_sub, 'fontsize', 12); %,'FontWeight','bold'
    correl=corr(max_A,my_var(var,:));
    title(strcat(kine_names(var),": ",string(round(correl,2))),'FontSize', 13)
end

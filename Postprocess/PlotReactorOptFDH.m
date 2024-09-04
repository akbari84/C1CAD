function PlotReactorOptFDH()
clc;

%% Parameters to set up path strings 
index=[1;2];
index_par=3;
n_index=size(index,1);
strReadFileMAT1=cell(n_index,1);
strReadFileMAT2=cell(n_index,1);

%% Path strings
strDataPath=Fun_DataPath();
strPlotPath=strcat(Fun_PlotPath(),'\ReactorOpt\rAcCoA');
for i=1:n_index
    strReadFileMAT1{i}=strcat(strDataPath,'\MAT\parameters-',num2str(index(i)),'.mat');
    strReadFileMAT2{i}=strcat(strDataPath,'\MAT\reactor_opt_par_',num2str(index_par),'-',num2str(index(i)),'.mat');
end
strPlotFilePDF1=strcat(strPlotPath,'\AssimRate_FDH_sp_par_',num2str(index_par),'.pdf');
strPlotFilePDF2=strcat(strPlotPath,'\kcat_FHD_sp_par_',num2str(index_par),'.pdf');
strPlotFilePDF3=strcat(strPlotPath,'\MW_FHD_sp_par_',num2str(index_par),'.pdf');

%% Memory allocation
sol_kcat=cell(n_index,1);
sol_rxn=cell(n_index,1);
sol_mc=cell(n_index,1);
sol_sc=cell(n_index,1);
assimilation_rate=zeros(n_index,2);
kcat=zeros(n_index,1);
MW=zeros(n_index,1);

%% Load data
for i=1:n_index
    load(strReadFileMAT1{i},'rxn');
    load(strReadFileMAT2{i},'multi_comp','single_comp','kapp');
    sol_kcat{i}=kapp;
    sol_rxn{i}=rxn;
    sol_mc{i}=multi_comp;
    sol_sc{i}=single_comp;
    clear('rxn','multi_comp','single_comp','kapp');
end

%% Parameters
j_fdh2=3;

%% Construct data matrices
for i=1:n_index
    assimilation_rate(i,1)=sol_mc{i}.assimilation_rate;
    assimilation_rate(i,2)=sol_sc{i}.assimilation_rate;
    kcat(i)=sol_kcat{i}(j_fdh2);
    MW(i)=sol_rxn{i}{j_fdh2}.MW;
end

%% Plot specifications 
x_label={'FDH';'Mo-FDH'};

BaseValue_As=1e-4;
BaseValue_kcat=1e-5;
BaseValue_MW=1e4;
color_mc_As=[140, 140, 140]/255;
color_sc_As=[200, 200, 200]/255;
color_kcat=[101, 200, 208]/255;
color_MW=[131, 139, 197]/255;

%% Plots
hh1=figure;
bar1=bar(assimilation_rate);
ylim([BaseValue_As,5e-2]);
set(gca,'XTick',1:n_index,'XTickLabel',x_label);
xtickangle(gca,90);
ylabel('$A_{\mathrm{CO}_2}$ (mol/s)','FontSize',20,'FontName','Arial','Interpreter','latex');
set(bar1(1),'BaseValue',BaseValue_As,'FaceColor',color_mc_As,'EdgeColor','black');
set(bar1(2),'BaseValue',BaseValue_As,'FaceColor',color_sc_As,'EdgeColor','black');
set(gca,'YScale','log','FontSize',20,'FontName','Arial');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperSize',[8 10]);
set(gcf,'PaperPosition',[0 0 8 10]);

hh2=figure;
bar2=bar(kcat);
ylim([BaseValue_kcat,1e2]);
set(gca,'XTick',1:n_index,'XTickLabel',x_label);
xtickangle(gca,90);
yticks([1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1 1e2]);
ylabel('$k_{\mathrm{cat}}$ (1/s)','FontSize',20,'FontName','Arial','Interpreter','latex');
set(bar2,'BaseValue',BaseValue_kcat,'FaceColor',color_kcat,'EdgeColor','black');
set(gca,'YScale','log','FontSize',20,'FontName','Arial');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperSize',[8 10]);
set(gcf,'PaperPosition',[0 0 8 10]);

hh3=figure;
bar3=bar(MW);
ylim([BaseValue_MW,4e5]);
set(gca,'XTick',1:n_index,'XTickLabel',x_label);
xtickangle(gca,90);
ylabel('$\mathrm{MW}$ (g/mol)','FontSize',20,'FontName','Arial','Interpreter','latex');
set(bar3,'BaseValue',BaseValue_MW,'FaceColor',color_MW,'EdgeColor','black');
set(gca,'YScale','log','FontSize',20,'FontName','Arial');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperSize',[8 10]);
set(gcf,'PaperPosition',[0 0 8 10]);

%% Print plots
print(hh1,'-dpdf',strPlotFilePDF1);
print(hh2,'-dpdf',strPlotFilePDF2);
print(hh3,'-dpdf',strPlotFilePDF3);
return
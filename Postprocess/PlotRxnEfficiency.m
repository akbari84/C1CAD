function PlotRxnEfficiency()
clc;

%% Parameters to set up path strings 
index=102;
index_par=1;

%% Path strings
strDataPath=Fun_DataPath();
if log10(index)<1
    strPlotPath=strcat(Fun_PlotPath(),'\ReactorOpt\rAcCoA');
elseif log10(index)>2
    strPlotPath=strcat(Fun_PlotPath(),'\ReactorOpt\rGly');
else
    return;
end
strReadFileMAT1=strcat(strDataPath,'\MAT\parameters-',num2str(index),'.mat');
strReadFileMAT2=strcat(strDataPath,'\MAT\reactor_opt_par_',num2str(index_par),'-',num2str(index),'.mat');
strPlotFilePDF1=strcat(strPlotPath,'\Therm_sp_par_',num2str(index_par),'-',num2str(index),'.pdf');
strPlotFilePDF2=strcat(strPlotPath,'\Therm_sp_km_par_',num2str(index_par),'-',num2str(index),'.pdf');
strPlotFilePDF3=strcat(strPlotPath,'\EnzSat_sp_km _par_',num2str(index_par),'-',num2str(index),'.pdf');

%% Load data
load(strReadFileMAT1,'n_rxn','rxn_id');
load(strReadFileMAT2,'multi_comp','single_comp','Xi','Xi1','Xi_','Xi1_','zeta_','zeta1_');

%% Parameters
kappa_atps_sat=12;
j_atps=n_rxn;

%% Memory allocation
gm=zeros(n_rxn,2);             % Column 1-> multi-compartment,    Column 2-> single-compartment
gm_=zeros(n_rxn,2);            % Column 1-> multi-compartment,    Column 2-> single-compartment
zt_=zeros(n_rxn,2);            % Column 1-> multi-compartment,    Column 2-> single-compartment

%% Calculate thermodynamic efficiencies
for j=1:n_rxn
    gm(j,1)=(1-Xi(j));
    gm(j,2)=(1-Xi1(j));
    gm_(j,1)=(1-Xi_(j));
    gm_(j,2)=(1-Xi1_(j));
    zt_(j,1)=zeta_(j);
    zt_(j,2)=zeta1_(j);
end
zt_(j_atps,1)=multi_comp.kappa_atps_/kappa_atps_sat;
zt_(j_atps,2)=single_comp.kappa_atps_/kappa_atps_sat;

%% Plot specifications 
gamma_BaseValue=0;
gamma_color_left=[45, 153, 181]/255;
gamma_color_left1=[140, 203, 219]/255;
gamma_color_left_=[45, 153, 181]/255;
gamma_color_left1_=[140, 203, 219]/255;

zeta_BaseValue=0;
zeta_color_left_=[71, 71, 71]/255;
zeta_color_left1_=[150, 150, 150]/255;

%% Plots
hh1=figure;
mu_factor=1e0;
bar1=bar(gm*mu_factor);
ylim([mu_factor*gamma_BaseValue,1]);
set(gca,'XTick',1:n_rxn,'XTickLabel',rxn_id);
xtickangle(gca,90);
ylabel('$\gamma$','FontSize',20,'FontName','Arial','Interpreter','latex');
set(bar1(1),'BaseValue',mu_factor*gamma_BaseValue,'FaceColor',gamma_color_left,'EdgeColor','black');
set(bar1(2),'BaseValue',mu_factor*gamma_BaseValue,'FaceColor',gamma_color_left1,'EdgeColor','black');
set(gca,'YColor',[0,0,0]/255,'FontSize',20,'FontName','Arial');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperSize',[20 12]);
set(gcf,'PaperPosition',[0 0 20 12]);

hh2=figure;
mu_factor=1e0;
bar1=bar(gm_*mu_factor);
ylim([mu_factor*gamma_BaseValue,1]);
set(gca,'XTick',1:n_rxn,'XTickLabel',rxn_id);
xtickangle(gca,90);
ylabel('$\gamma$','FontSize',20,'FontName','Arial','Interpreter','latex');
set(bar1(1),'BaseValue',mu_factor*gamma_BaseValue,'FaceColor',gamma_color_left_,'EdgeColor','black');
set(bar1(2),'BaseValue',mu_factor*gamma_BaseValue,'FaceColor',gamma_color_left1_,'EdgeColor','black');
set(gca,'YColor',[0,0,0]/255,'FontSize',20,'FontName','Arial');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperSize',[20 12]);
set(gcf,'PaperPosition',[0 0 20 12]);

hh3=figure;
mu_factor=1e0;
bar1=bar(zt_*mu_factor);
ylim([mu_factor*zeta_BaseValue,1]);
set(gca,'XTick',1:n_rxn,'XTickLabel',rxn_id);
xtickangle(gca,90);
ylabel('$\zeta$','FontSize',20,'FontName','Arial','Interpreter','latex');
set(bar1(1),'BaseValue',mu_factor*zeta_BaseValue,'FaceColor',zeta_color_left_,'EdgeColor','black');
set(bar1(2),'BaseValue',mu_factor*zeta_BaseValue,'FaceColor',zeta_color_left1_,'EdgeColor','black');
set(gca,'YColor',[0,0,0]/255,'FontSize',20,'FontName','Arial');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperSize',[20 12]);
set(gcf,'PaperPosition',[0 0 20 12]);

%% Print plots
print(hh1,'-dpdf',strPlotFilePDF1);
print(hh2,'-dpdf',strPlotFilePDF2);
print(hh3,'-dpdf',strPlotFilePDF3);

return



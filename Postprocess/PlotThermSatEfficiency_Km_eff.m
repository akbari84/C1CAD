function PlotThermSatEfficiency_Km_eff()
clc;

%% Parameters to set up path strings 
index=102;
index_par=1;

%% Path strings
strDataPath=Fun_DataPath();
if log10(index)<1
    strPlotPath=strcat(Fun_PlotPath(),'\ReactorOpt\rAcCoA');
end
if log10(index)>2
    strPlotPath=strcat(Fun_PlotPath(),'\ReactorOpt\rGly');
end
strReadFileMAT1=strcat(strDataPath,'\MAT\parameters-',num2str(index),'.mat');
strReadFileMAT2=strcat(strDataPath,'\MAT\reactor_opt_par_',num2str(index_par),'-',num2str(index),'.mat');
strPlotFilePDF1=strcat(strPlotPath,'\ThermSatEffic_mc_km_eff_par_',num2str(index_par),'-',num2str(index),'.pdf');
strPlotFilePDF2=strcat(strPlotPath,'\ThermSatEffic_sc_km_eff_par_',num2str(index_par),'-',num2str(index),'.pdf');
strPlotFilePDF3=strcat(strPlotPath,'\AssimRate_km_eff_par_',num2str(index_par),'-',num2str(index),'.pdf');

%% Load data
load(strReadFileMAT1,'rxn','n_met','n_rxn','S','Km','rxn_id');
load(strReadFileMAT2,'multi_comp','single_comp','Xi','Xi1','Xi_','Xi1_','kapp','kapp1');

%% Parameters
if log10(index)<1
    jc=[1;2;3;4;5;6;7;8;9;10];
end
if log10(index)>2
    jc=[1;2;3;4;5;6;7;8;9;10];
end
n_rxn_c=size(jc,1);
ar_photosynthesis=0.00044718;           % mol/s for N. gaditana

%% Memory allocation
gamma=zeros(n_rxn_c,1);
gamma1=zeros(n_rxn_c,1);
gamma_=zeros(n_rxn_c,1);
gamma1_=zeros(n_rxn_c,1);
zeta=zeros(n_rxn_c,1);
zeta1=zeros(n_rxn_c,1);
zeta_=zeros(n_rxn_c,1);
zeta1_=zeros(n_rxn_c,1);

%% Unit conversion for Km
Km=Km*1e-6;         % Micro mole in Excel sheet -> mol in this routine 

%% Correct for reaction directions
for j=1:n_rxn
    S(:,j)=S(:,j)*rxn{j}.dir;
end

%% Turnover number
kcat=kapp(jc);
kcat1=kapp1(jc);

%% Calculate thermodynamic efficiency
for j=1:n_rxn_c
    gamma(j)=1-Xi(jc(j));
    gamma1(j)=1-Xi1(jc(j));
    gamma_(j)=1-Xi_(jc(j));
    gamma1_(j)=1-Xi1_(jc(j));
end

%% Calculate saturation efficiency
for j=1:n_rxn_c
    zeta(j)=Fun_Zeta(multi_comp.xc,Km,S,n_met,jc(j));
    zeta1(j)=Fun_Zeta(single_comp.xc,Km,S,n_met,jc(j));
    zeta_(j)=Fun_Zeta(multi_comp.xc,Km,S,n_met,jc(j));
    zeta1_(j)=Fun_Zeta(single_comp.xc,Km,S,n_met,jc(j));
end

%% Plot specifications 
eff_km_color=[33, 33, 33]/255;
eff_km_color_s=[150, 150, 150]/255;         % shadow color
eff_km_text_color=[0, 0, 0]/255;
eff_photosynthesis_color=[0, 196, 32]/255;
eff_color=[255, 200, 3]/255;
eff_color_s=[255, 248, 222]/255;              % shadow color
eff_text_color=[138, 79, 12]/255;

eff_edge_color=[0,0,0]/255;
eff_edge_color_s=[50,50,50]/255;

LW=1;
sz=36;
rr_s=0.01;                              % radius of shadow
data_label=rxn_id(jc);
dd=0.005;

BaseValue_Ar=1e-4;
%BaseValue_Ar=1e0;
assim_rate=[multi_comp.assimilation_rate_ multi_comp.assimilation_rate-multi_comp.assimilation_rate_;single_comp.assimilation_rate_ single_comp.assimilation_rate-single_comp.assimilation_rate_];
ar_color_bottom=[33, 33, 33]/255;
ar_color_top=[255, 200, 3]/255;
x_label={'MC';'SC'};

%% Plots
hh1=figure;
scatter3(gamma,zeta,kcat,sz,'o','MarkerFaceColor',eff_color,'MarkerEdgeColor',eff_edge_color,'LineWidth',LW);
hold on
for j=1:n_rxn_c
    rectangle('Position',[gamma(j)-rr_s zeta(j)-rr_s rr_s rr_s],'Curvature',1,'FaceColor',eff_color_s,'EdgeColor',eff_edge_color_s,'LineWidth',LW);
    hold on
end

scatter3(gamma_,zeta_,kcat,sz,'o','MarkerFaceColor',eff_km_color,'MarkerEdgeColor',eff_edge_color,'LineWidth',LW);
hold on
for j=1:n_rxn_c
    rectangle('Position',[gamma_(j)-rr_s zeta_(j)-rr_s rr_s rr_s],'Curvature',1,'FaceColor',eff_km_color_s,'EdgeColor',eff_edge_color_s,'LineWidth',LW);
    hold on
end

for j=1:n_rxn_c
    text(gamma(j)+dd,zeta(j),kcat(j),data_label{j},'Color',eff_text_color);
end
for j=1:n_rxn_c
    text(gamma_(j)+dd,zeta_(j),kcat(j),data_label{j},'Color',eff_km_text_color);
end

xlabel('Thermodynamic efficiency','FontSize',20,'FontName','Arial','Interpreter','latex');
ylabel('Saturation efficiency','FontSize',20,'FontName','Arial','Interpreter','latex');
set(gca,'ZScale','log','FontSize',20,'FontName','Arial');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperSize',[10 10]);
set(gcf,'PaperPosition',[0 0 10 10]);

hh2=figure;
scatter(gamma1,zeta1,sz,'o','MarkerFaceColor',eff_color,'MarkerEdgeColor',eff_edge_color,'LineWidth',LW);
for j=1:n_rxn_c
    text(gamma1(j)+dd,zeta1(j),data_label{j},'Color',eff_text_color);
end
hold on
scatter(gamma1_,zeta1_,sz,'o','MarkerFaceColor',eff_km_color,'MarkerEdgeColor',eff_edge_color,'LineWidth',LW);
for j=1:n_rxn_c
    text(gamma1_(j)+dd,zeta1_(j),data_label{j},'Color',eff_km_text_color);
end

xlabel('Thermodynamic efficiency','FontSize',20,'FontName','Arial','Interpreter','latex');
ylabel('Saturation efficiency','FontSize',20,'FontName','Arial','Interpreter','latex');
set(gca,'FontSize',20,'FontName','Arial');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperSize',[10 10]);
set(gcf,'PaperPosition',[0 0 10 10]);

hh3=figure;
bar1=bar(assim_rate,'stacked');
hold on 
yline(ar_photosynthesis,'Color',eff_photosynthesis_color);
%ylim([BaseValue_Ar,5e-2]);
ylim([BaseValue_Ar,1e1]);
set(gca,'XTick',1:2,'XTickLabel',x_label);
xtickangle(gca,90);
ylabel('$A_{\mathrm{CO}_2}$ (mol/s)','FontSize',20,'FontName','Arial','Interpreter','latex');
set(bar1(1),'BaseValue',BaseValue_Ar,'FaceColor',ar_color_bottom,'EdgeColor','black');
set(bar1(2),'BaseValue',BaseValue_Ar,'FaceColor',ar_color_top,'EdgeColor','black');
set(gca,'YScale','log','FontSize',20,'FontName','Arial');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperSize',[8 10]);
set(gcf,'PaperPosition',[0 0 8 10]);

%% Print plots
print(hh1,'-dpdf',strPlotFilePDF1);
print(hh2,'-dpdf',strPlotFilePDF2);
print(hh3,'-dpdf',strPlotFilePDF3);
return

function zeta=Fun_Zeta(x,Km,S,n_met,j)
prod_s=1;
prod_p=1;
n_s=0;
n_p=0;
for i=1:n_met
    if S(i,j)==0 || Km(i,j)==0
        continue;
    elseif S(i,j)>0
        prod_p=prod_p*(x(i)/Km(i,j))^S(i,j);
        n_p=n_p+1;
    else
        prod_s=prod_s*(x(i)/Km(i,j))^(-S(i,j));
        n_s=n_s+1;
    end
end
if n_s==0
    prod_s=0;
end
if n_p==0
    prod_p=0;
end
zeta=prod_s/(1+prod_s+prod_p);
return
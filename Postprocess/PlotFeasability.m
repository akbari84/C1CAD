function PlotFeasability()
clc;

%% Parameters to set up path strings 
index_par=1;

%% Path strings
strPlotPath=strcat(Fun_PlotPath(),'\ReactorOpt_swp');
strPlotFilePDF1=strcat(strPlotPath,'\rAcCoA\FeasRegion_rAcCoA_sc_ub_par_',num2str(index_par),'.pdf');
strPlotFilePDF2=strcat(strPlotPath,'\rAcCoA\FeasRegion_rAcCoA_mc_ub_par_',num2str(index_par),'.pdf');
strPlotFilePDF3=strcat(strPlotPath,'\rAcCoA\FeasRegion_rAcCoA_sc_H2_par_',num2str(index_par),'.pdf');
strPlotFilePDF4=strcat(strPlotPath,'\rAcCoA\FeasRegion_rAcCoA_mc_H2_par_',num2str(index_par),'.pdf');
strPlotFilePDF5=strcat(strPlotPath,'\rGly\FeasRegion_rGly_sc_ub_par_',num2str(index_par),'.pdf');
strPlotFilePDF6=strcat(strPlotPath,'\rGly\FeasRegion_rGly_mc_ub_par_',num2str(index_par),'.pdf');
strPlotFilePDF7=strcat(strPlotPath,'\rGly\FeasRegion_rGly_sc_H2_par_',num2str(index_par),'.pdf');
strPlotFilePDF8=strcat(strPlotPath,'\rGly\FeasRegion_rGly_mc_H2_par_',num2str(index_par),'.pdf');
strPlotFilePDF9=strcat(strPlotPath,'\rAcCoA\uc_rAcCoA_sc_H2_par_',num2str(index_par),'.pdf');
strPlotFilePDF10=strcat(strPlotPath,'\rAcCoA\uc_rAcCoA_mc_H2_par_',num2str(index_par),'.pdf');
strPlotFilePDF11=strcat(strPlotPath,'\rGly\uc_rGly_sc_H2_par_',num2str(index_par),'.pdf');
strPlotFilePDF12=strcat(strPlotPath,'\rGly\uc_rGly_mc_H2_par_',num2str(index_par),'.pdf');

%% Parameters
np_pHc=81;
pHc_min=5.5;
pHc_max=8.5;

np_H2=81;
x_c_H2_min=0.02;
x_c_H2_max=0.1;

np_ub=81;
x_ub_min=0.05;
x_ub_max=0.15;

x_c_H2_0=0.05;
x_ub_0=0.05;
beta0=0.4;

%% Memory allocation
is_feas_raccoa_sc_ub=zeros(np_pHc,np_ub);
is_feas_raccoa_mc_ub=zeros(np_pHc,np_ub);
is_feas_raccoa_sc_H2=zeros(np_pHc,np_H2);
is_feas_raccoa_mc_H2=zeros(np_pHc,np_H2);
is_feas_rgly_sc_ub=zeros(np_pHc,np_ub);
is_feas_rgly_mc_ub=zeros(np_pHc,np_ub);
is_feas_rgly_sc_H2=zeros(np_pHc,np_H2);
is_feas_rgly_mc_H2=zeros(np_pHc,np_H2);
uc_raccoa_sc_H2=zeros(np_pHc,np_H2);
uc_raccoa_mc_H2=zeros(np_pHc,np_H2);
uc_rgly_sc_H2=zeros(np_pHc,np_H2);
uc_rgly_mc_H2=zeros(np_pHc,np_H2);

%% Grid generation
pHc=linspace(pHc_min,pHc_max,np_pHc);
x_c_H2=linspace(x_c_H2_min,x_c_H2_max,np_H2);
x_ub=linspace(x_ub_min,x_ub_max,np_ub);

%% Calculate feasibility
index=2;
for i=1:np_pHc
    for j=1:np_ub
        [is_feas_raccoa_sc_ub(i,j),~,~]=Fun_FeasabilitySC_rAcCoA(pHc(i),beta0,x_c_H2_0,x_ub(j),index);
    end
    for j=1:np_H2
        [is_feas_raccoa_sc_H2(i,j),uc_raccoa_sc_H2(i,j),~]=Fun_FeasabilitySC_rAcCoA(pHc(i),beta0,x_c_H2(j),x_ub_0,index);
    end
end
for i=1:np_pHc
    for j=1:np_ub
        [is_feas_raccoa_mc_ub(i,j),~]=Fun_FeasabilityMC_rAcCoA(pHc(i),beta0,x_c_H2_0,x_ub(j),index);
    end
    for j=1:np_H2
        [is_feas_raccoa_mc_H2(i,j),uc_raccoa_mc_H2(i,j)]=Fun_FeasabilityMC_rAcCoA(pHc(i),beta0,x_c_H2(j),x_ub_0,index);
    end
end

index=102;
for i=1:np_pHc
    for j=1:np_ub
        [is_feas_rgly_sc_ub(i,j),~,~]=Fun_FeasabilitySC_rGly(pHc(i),beta0,x_c_H2_0,x_ub(j),index);
    end
    for j=1:np_H2
        [is_feas_rgly_sc_H2(i,j),uc_rgly_sc_H2(i,j),~]=Fun_FeasabilitySC_rGly(pHc(i),beta0,x_c_H2(j),x_ub_0,index);
    end
end
for i=1:np_pHc
    for j=1:np_ub
        [is_feas_rgly_mc_ub(i,j),~]=Fun_FeasabilityMC_rGly(pHc(i),beta0,x_c_H2_0,x_ub(j),index);
    end
    for j=1:np_H2
        [is_feas_rgly_mc_H2(i,j),uc_rgly_mc_H2(i,j)]=Fun_FeasabilityMC_rGly(pHc(i),beta0,x_c_H2(j),x_ub_0,index);
    end
end

%% Plots
hh1=figure;
[X,Y]=meshgrid(pHc,x_ub);
contourf(X,Y,is_feas_raccoa_sc_ub',[1 1]);
xlabel('$\mathrm{pH}_{\mathrm{c}}$','FontSize',20,'FontName','Times New Roman','Interpreter','latex');
ylabel('$x^{\mathrm{ub}}$ (M)','FontSize',20,'FontName','Times New Roman','Interpreter','latex');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperSize',[10 10]);
set(gcf,'PaperPosition',[0 0 10 10]);

hh2=figure;
contourf(X,Y,is_feas_raccoa_mc_ub',[1 1]);
xlabel('$\mathrm{pH}_{\mathrm{c}}$','FontSize',20,'FontName','Times New Roman','Interpreter','latex');
ylabel('$x^{\mathrm{ub}}$ (M)','FontSize',20,'FontName','Times New Roman','Interpreter','latex');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperSize',[10 10]);
set(gcf,'PaperPosition',[0 0 10 10]);

hh3=figure;
[X,Y]=meshgrid(pHc,x_c_H2);
contourf(X,Y,is_feas_raccoa_sc_H2',[1 1]);
xlabel('$\mathrm{pH}_{\mathrm{c}}$','FontSize',20,'FontName','Times New Roman','Interpreter','latex');
ylabel('$x^c_{\mathrm{H}_2}$ (M)','FontSize',20,'FontName','Times New Roman','Interpreter','latex');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperSize',[10 10]);
set(gcf,'PaperPosition',[0 0 10 10]);

hh4=figure;
contourf(X,Y,is_feas_raccoa_mc_H2',[1 1]);
xlabel('$\mathrm{pH}_{\mathrm{c}}$','FontSize',20,'FontName','Times New Roman','Interpreter','latex');
ylabel('$x^c_{\mathrm{H}_2}$ (M)','FontSize',20,'FontName','Times New Roman','Interpreter','latex');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperSize',[10 10]);
set(gcf,'PaperPosition',[0 0 10 10]);

hh5=figure;
[X,Y]=meshgrid(pHc,x_ub);
contourf(X,Y,is_feas_rgly_sc_ub',[1 1]);
xlabel('$\mathrm{pH}_{\mathrm{c}}$','FontSize',20,'FontName','Times New Roman','Interpreter','latex');
ylabel('$x^{\mathrm{ub}}$ (M)','FontSize',20,'FontName','Times New Roman','Interpreter','latex');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperSize',[10 10]);
set(gcf,'PaperPosition',[0 0 10 10]);

hh6=figure;
contourf(X,Y,is_feas_rgly_mc_ub',[1 1]);
xlabel('$\mathrm{pH}_{\mathrm{c}}$','FontSize',20,'FontName','Times New Roman','Interpreter','latex');
ylabel('$x^{\mathrm{ub}}$ (M)','FontSize',20,'FontName','Times New Roman','Interpreter','latex');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperSize',[10 10]);
set(gcf,'PaperPosition',[0 0 10 10]);

hh7=figure;
[X,Y]=meshgrid(pHc,x_c_H2);
contourf(X,Y,is_feas_rgly_sc_H2',[1 1]);
xlabel('$\mathrm{pH}_{\mathrm{c}}$','FontSize',20,'FontName','Times New Roman','Interpreter','latex');
ylabel('$x^c_{\mathrm{H}_2}$ (M)','FontSize',20,'FontName','Times New Roman','Interpreter','latex');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperSize',[10 10]);
set(gcf,'PaperPosition',[0 0 10 10]);

hh8=figure;
contourf(X,Y,is_feas_rgly_mc_H2',[1 1]);
xlabel('$\mathrm{pH}_{\mathrm{c}}$','FontSize',20,'FontName','Times New Roman','Interpreter','latex');
ylabel('$x^c_{\mathrm{H}_2}$ (M)','FontSize',20,'FontName','Times New Roman','Interpreter','latex');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperSize',[10 10]);
set(gcf,'PaperPosition',[0 0 10 10]);

hh9=figure;
[X,Y]=meshgrid(pHc,x_c_H2);
contourf(X,Y,uc_raccoa_sc_H2');
colorbar;
xlabel('$\mathrm{pH}_{\mathrm{c}}$','FontSize',20,'FontName','Times New Roman','Interpreter','latex');
ylabel('$x^c_{\mathrm{H}_2}$ (M)','FontSize',20,'FontName','Times New Roman','Interpreter','latex');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperSize',[10 10]);
set(gcf,'PaperPosition',[0 0 10 10]);

hh10=figure;
contourf(X,Y,uc_raccoa_mc_H2');
colorbar;
xlabel('$\mathrm{pH}_{\mathrm{c}}$','FontSize',20,'FontName','Times New Roman','Interpreter','latex');
ylabel('$x^c_{\mathrm{H}_2}$ (M)','FontSize',20,'FontName','Times New Roman','Interpreter','latex');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperSize',[10 10]);
set(gcf,'PaperPosition',[0 0 10 10]);

hh11=figure;
[X,Y]=meshgrid(pHc,x_c_H2);
contourf(X,Y,uc_rgly_sc_H2');
colorbar;
xlabel('$\mathrm{pH}_{\mathrm{c}}$','FontSize',20,'FontName','Times New Roman','Interpreter','latex');
ylabel('$x^c_{\mathrm{H}_2}$ (M)','FontSize',20,'FontName','Times New Roman','Interpreter','latex');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperSize',[10 10]);
set(gcf,'PaperPosition',[0 0 10 10]);

hh12=figure;
contourf(X,Y,uc_rgly_mc_H2');
colorbar;
xlabel('$\mathrm{pH}_{\mathrm{c}}$','FontSize',20,'FontName','Times New Roman','Interpreter','latex');
ylabel('$x^c_{\mathrm{H}_2}$ (M)','FontSize',20,'FontName','Times New Roman','Interpreter','latex');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperSize',[10 10]);
set(gcf,'PaperPosition',[0 0 10 10]);

%% Print plots
print(hh1,'-dpdf',strPlotFilePDF1);
print(hh2,'-dpdf',strPlotFilePDF2);
print(hh3,'-dpdf',strPlotFilePDF3);
print(hh4,'-dpdf',strPlotFilePDF4);
print(hh5,'-dpdf',strPlotFilePDF5);
print(hh6,'-dpdf',strPlotFilePDF6);
print(hh7,'-dpdf',strPlotFilePDF7);
print(hh8,'-dpdf',strPlotFilePDF8);
print(hh9,'-dpdf',strPlotFilePDF9);
print(hh10,'-dpdf',strPlotFilePDF10);
print(hh11,'-dpdf',strPlotFilePDF11);
print(hh12,'-dpdf',strPlotFilePDF12);
return

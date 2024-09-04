function ReactorOpt_v3()
% Reductive Glycine
clc;

%% Parameters to set up path strings 
index=102;
index_par=3;       

%% Path strings
strDataPath=Fun_DataPath();
strPlotPath=strcat(Fun_PlotPath(),'\ReactorOpt\rGly');
strReadFileMAT1=strcat(strDataPath,'\MAT\parameters-',num2str(index),'.mat');
strPlotFilePDF1=strcat(strPlotPath,'\EnzymeConc_sp_par_',num2str(index_par),'-',num2str(index),'.pdf');
strPlotFilePDF2=strcat(strPlotPath,'\Flux_sp_par_',num2str(index_par),'-',num2str(index),'.pdf');
strPlotFilePDF3=strcat(strPlotPath,'\Conc_qc_sp_par_',num2str(index_par),'-',num2str(index),'.pdf');
strPlotFilePDF4=strcat(strPlotPath,'\EnzymeConc_sp_km_par_',num2str(index_par),'-',num2str(index),'.pdf');
strPlotFilePDF5=strcat(strPlotPath,'\Flux_sp_km_par_',num2str(index_par),'-',num2str(index),'.pdf');
strPlotFilePDF6=strcat(strPlotPath,'\Conc_qc_sp_km_par_',num2str(index_par),'-',num2str(index),'.pdf');  
strPlotFilePDF7=strcat(strPlotPath,'\EnzymeConc_sp_mc_km_eff_par_',num2str(index_par),'-',num2str(index),'.pdf');
strPlotFilePDF8=strcat(strPlotPath,'\EnzymeConc_sp_sc_km_eff_par_',num2str(index_par),'-',num2str(index),'.pdf');
strPlotFilePDF9=strcat(strPlotPath,'\Flux_sp_mc_km_eff_par_',num2str(index_par),'-',num2str(index),'.pdf');
strPlotFilePDF10=strcat(strPlotPath,'\Flux_sp_sc_km_eff_par_',num2str(index_par),'-',num2str(index),'.pdf');
strPlotFilePDF11=strcat(strPlotPath,'\AssimRate_km_eff_par_',num2str(index_par),'-',num2str(index),'.pdf');
strSaveFileMAT=strcat(strDataPath,'\MAT\reactor_opt_par_',num2str(index_par),'-',num2str(index),'.mat');

%% Load data
load(strReadFileMAT1,'met','rxn','n_met','n_rxn','S','Km','met_id','rxn_id');

%% Constants
Rgas=0.00831446261815324;                           % kJ/K/mol   
F=96.48533212;                                      % kC/mol (or kJ/mol/V)
A_PHI_PZ=0.392139132830784;  % Alberty at 25 C      % mol^(-1/2) kg^(1/2)
B_PHI_PZ=1.2;                                       % mol^(-1/2) kg^(1/2)

%% Indices
i_h2=1;
i_Fdrd=3;
i_Fdox=4;
i_TRXrd=5;
i_TRXox=6;
i_co2=2;
i_pyr=27;
i_nadp=8;
i_nadph=9;
i_nad=10;
i_nadh=11;
i_atp=12;
i_adp=13;
i_amp=14;
i_pi=15;
i_ppi=16;

i_rd=[i_nadh;i_nadph;i_Fdrd;i_TRXrd];
i_ox=[i_nad;i_nadp;i_Fdox;i_TRXox];

comp_j=[2;2;2;2;2;2;2;2;2;2;3;3;3;3;1;1;1];
jc=[1;2;3;4;5;6;7;8;9;10];
jc_f=[1;2;3;4;5;6;7;8;9];
jc_m=10;
jf=[11;12;13;14];
je_v=[15;16];
je_a=17;
j_fdh1=1;
j_fdh2=2;
j_atps=17;
j_nadh=11;              % HOX
j_nadph=12;             % HND
j_fdrd=13;              % HYD 
j_trxrd=14;             % FRH 
n_rxn_c=size(jc,1);
n_rxn_f=size(jf,1);
n_rxn_ev=size(je_v,1);
n_rxn1=n_rxn_c+n_rxn_f+n_rxn_ev;
n_pH=5;                 % 3 compartments for multicompartment case + 1 for single compartment case + 1 reference 

%% Memory allocation
DfG0tr=zeros(n_met,n_pH);
DrG0tr=zeros(n_rxn,1);
DrG0tr1=zeros(n_rxn,1);
DrG0tr_pH7=zeros(n_rxn,1);
Keq=zeros(n_rxn,1);
Keq1=zeros(n_rxn,1);
DrG_H=zeros(n_pH,1);
DrG0tr_S=zeros(n_rxn,n_pH);
DfG0tr_w=zeros(n_pH,1);
f_gamma=zeros(n_pH,1);
Ln_gamma_H=zeros(n_pH,1);
gamma_H=zeros(n_pH,1);
Ln_gamma_w=zeros(n_pH,1);
x_CO2_Q=zeros(n_met,1);
x_H2_f=zeros(n_met,1);
varth=zeros(n_rxn,1);

xq=zeros(n_met,1); 
%xc=zeros(n_met,1); 
xF=zeros(n_met,1); 
xE=zeros(n_met,1); 
xf=zeros(n_met,1); 
xe=zeros(n_met,1); 
xC=zeros(n_met,1); 
Enz=zeros(n_rxn,1);
v=zeros(n_rxn,1);
kapp=zeros(n_rxn,1);
w=zeros(n_rxn_c,1);
Xi=zeros(n_rxn,1);

xq_=zeros(n_met,1); 
xc_=zeros(n_met,1); 
xF_=zeros(n_met,1); 
xE_=zeros(n_met,1); 
xf_=zeros(n_met,1); 
xe_=zeros(n_met,1); 
xC_=zeros(n_met,1); 
Enz_=zeros(n_rxn,1);
v_=zeros(n_rxn,1);
Xi_=zeros(n_rxn,1);
zeta_=zeros(n_rxn,1);

xq1=zeros(n_met,1); 
xc1=zeros(n_met,1); 
xR1=zeros(n_met,1); 
xC1=zeros(n_met,1); 
Enz1=zeros(n_rxn,1);
v1=zeros(n_rxn,1);
kapp1=zeros(n_rxn,1);
w1=zeros(n_rxn1,1);
Xi1=zeros(n_rxn,1);

xq1_=zeros(n_met,1); 
xc1_=zeros(n_met,1); 
xR1_=zeros(n_met,1); 
xC1_=zeros(n_met,1); 
Enz1_=zeros(n_rxn,1);
v1_=zeros(n_rxn,1);
Xi1_=zeros(n_rxn,1);
zeta1_=zeros(n_rxn,1);

%% Parameters
T=25+273.15;                                        % K
pH=[7;6;8.0;7.5;7];                                 % -             1->e, 2->c, 3->f, 4->single compartment,  5->reference for flux dir correction        
pH_p=2;                                             % -
Is=[0.1;0.1;0.1;0.1;0.1];                           % M             1->e, 2->c, 3->f, 4->single compartment,  5->reference for flux dir correction 
Is_p=0.1;                                           % M 

DPsi=-0.1;                                          % V 
DE_FD=-0.399;                                       % V             DE_FD:=E_FDrd-E_FDox     
DG00w=-238.7;                                       % kJ/mol 
Enz_tot_v=340;                                      % g/L           Based on E. coli cytoplasmatic density 
Enz_tot_a=2.6569e-09;                               % mol/m^2       Estimated based on the diameters of ATPS, assuming each protein cannot fit in a square smaller than 25x25 nm 

x_CO2_Q(i_co2)=0.05;                                % M           The concentration of species in the CO2 feed stream (i.e. Q); solubality at 25C = 0.003  
x_H2_f(i_h2)=0.05;                                  % M           The concentration of species in the H2 feed stream (i.e. H2 to f-reactor); solubality at 25C = 0.001 
x_c_H2=[5e-2,5e-2];                                 % M           The concetration of H2 in the c stream [x_c_H2(1)->multi-comp, x_c_H2(2)->single comp]. We can treat it as a free parameter without affecting mass-balances around c,f compartments (unreacted H2 circulates). 
V_c=1.0e3;                                          % L
V1=1.0e3;                                           % L
beta_f=0.2;                                         % -           beta_f:=FF/c  
beta_e=0.2;                                         % -           beta_e:=EE/c
beta_r=0.4;                                         % -           beta_r:=RR/c
eta_hat=1;                                          % -           Degree of separation for unit operation at outlet node: fl('FF')*x(i_x,'FF')+fl('EE')*x(i_x,'EE')=eta_hat*fl('c')*x(i_x,'c')
eta=0.99;                                           % -           Degree of separation for unit operation at outlet node: fl('CC')*x(pyr,'CC')=eta*fl('c')*x(pyr,'c')   

%% ATPS parameters
kon_pi=100;                                         % M^-1 s^-1
kon_adp=4.0e7;                                      % M^-1 s^-1
kon_atp=4.0e8;                                      % M^-1 s^-1
kon_hp=4.0e11;                                      % M^-1 s^-1
kon_hc=4.0e11;                                      % M^-1 s^-1
k_dt=1e3;                                           % s^-1

Kd_pi=2;                                            % M
Kd_adp=9.0e-6;                                      % M
Kd_atp=5.0e-3;                                      % M
Kd_hp1=2.0e-7;                                      % M
Kd_hc1=2.0e-7;                                      % M
Keq_dt=10;                                          % -

omega_p1=10;
omega_p2=10^0.2;
omega_c1=10;
omega_c2=0.1;

Kd_hp2=omega_p1*Kd_hp1;                             % M
Kd_hc2=omega_c1*Kd_hc1;                             % M
Kd_hp3=omega_p2*Kd_hp2;                             % M
Kd_hc3=omega_c2*Kd_hc2;                             % M
n_step=10;

atps_par=[kon_pi;kon_adp;kon_atp;kon_hp;kon_hc;k_dt;Kd_pi;Kd_adp;Kd_atp;Kd_hp1;Kd_hc1;Keq_dt;Kd_hp2;Kd_hc2;Kd_hp3;Kd_hc3;n_step];

%% Dependent parameters
DG00_FD=-F*DE_FD;                                   % kJ/mol  
x_H=10.^(-pH);                                      % M 
x_H_p=10^(-pH_p);                                   % M 

%% Cofactor proteins (parametrized)
met{i_Fdrd}.DfG00_A=DG00_FD;                        % Note: we always choose met{i_FDox}.DfG00_A=0

%% Dissociation constants
for i=1:n_met
    met{i}.K=10.^(-met{i}.pKa);
    met{i}.KK=ones(1,met{i}.N+1);
    for k=2:met{i}.N+1
        met{i}.KK(k)=met{i}.KK(k-1)*met{i}.K(k-1);
    end
end

%% Activity coefficients
for i=1:n_pH
    f_gamma(i)=f_gamma_EDH(Is(i),A_PHI_PZ,B_PHI_PZ);
    Ln_gamma_H(i)=-f_gamma(i);
    Ln_gamma_w(i)=0;                            % To exactly match Alberty's water activity
    gamma_H(i)=exp(Ln_gamma_H(i)); 
end
f_gamma_p=f_gamma_EDH(Is_p,A_PHI_PZ,B_PHI_PZ);
Ln_gamma_H_p=-f_gamma_p;
                                   
%% Formation energies
for i_c=1:n_pH
    for i=1:n_met
        sum1=0;
        for k=1:met{i}.N+1
            if met{i}.is_aq==1 
                zk=met{i}.z_A+k-1;
            else
                zk=0;
            end
            Ln_gamma_k=-zk^2*f_gamma(i_c);
            gamma_k=exp(Ln_gamma_k);
            sum1=sum1+gamma_H(i_c)^(k-1)/10^((k-1)*pH(i_c))/gamma_k/met{i}.KK(k);
        end 
        DfG0tr(i,i_c)=met{i}.DfG00_A+Rgas*T*(met{i}.NH_A*(log(10)*pH(i_c)-Ln_gamma_H(i_c))-log(sum1));
    end
end

%% Reaction energies
for i_c=1:n_pH
    DrG_H(i_c)=Rgas*T*(Ln_gamma_H(i_c)-Ln_gamma_H_p-log(10)*(pH(i_c)-pH_p))+F*DPsi;                         % Proton translocation contribution
    DrG0tr_S(:,i_c)=S'*DfG0tr(:,i_c);                                                                       % Stoichiometry contribution
    DfG0tr_w(i_c)=DG00w+Rgas*T*(Ln_gamma_w(i_c)-2*Ln_gamma_H(i_c)+2*log(10)*pH(i_c));                       % Water contribution
end

for j=1:n_rxn
    DrG0tr(j)=DrG0tr_S(j,comp_j(j))+rxn{j}.nH_transport*DrG_H(comp_j(j))+rxn{j}.S_w*DfG0tr_w(comp_j(j));
    Keq(j)=exp(-DrG0tr(j)/Rgas/T);
    DrG0tr1(j)=DrG0tr_S(j,4)+rxn{j}.nH_transport*DrG_H(4)+rxn{j}.S_w*DfG0tr_w(4);
    DrG0tr_pH7(j)=DrG0tr_S(j,5)+rxn{j}.nH_transport*DrG_H(5)+rxn{j}.S_w*DfG0tr_w(5);
    Keq1(j)=exp(-DrG0tr1(j)/Rgas/T);
end

%% Unit conversion for Km
Km=Km*1e-6;         % Micro mole in Excel sheet -> mol in this routine 

%% Correct for reaction directions
for j=1:n_rxn-1
    prd1=1;
    for i=1:n_met
        if Km(i,j)~=0 
            prd1=prd1*Km(i,j)^S(i,j);
        end
    end
    if rxn{j}.dir==1 && rxn{j}.dir_kin==-1
        rxn{j}.kcat=rxn{j}.kcat*exp(-DrG0tr_pH7(j)/Rgas/T)/prd1;
    end
    if rxn{j}.dir==-1 && rxn{j}.dir_kin==1
        rxn{j}.kcat=rxn{j}.kcat*exp(DrG0tr_pH7(j)/Rgas/T)*prd1;
    end
end
for j=1:n_rxn
    DrG0tr(j)=DrG0tr(j)*rxn{j}.dir;
    DrG0tr1(j)=DrG0tr1(j)*rxn{j}.dir;
    Keq(j)=Keq(j)^rxn{j}.dir;
    Keq1(j)=Keq1(j)^rxn{j}.dir;
    S(:,j)=S(:,j)*rxn{j}.dir;
    rxn{j}.S_w=rxn{j}.S_w*rxn{j}.dir;
    rxn{j}.nH_transport=rxn{j}.nH_transport*rxn{j}.dir;
end

%% Apparent rate constants
for j=1:size(jc_f,1)
    kapp(jc_f(j))=rxn{jc_f(j)}.kcat;
end
for j=1:size(jc_m,1)
    kapp(jc_m(j))=rxn{jc_m(j)}.kcat;
end
for j=1:size(jf,1)
    kapp(jf(j))=rxn{jf(j)}.kcat;
end
for j=1:size(je_v,1)
    kapp(je_v(j))=rxn{je_v(j)}.kcat;
end
for j=1:size(je_a,1)
    kapp(je_a(j))=rxn{je_a(j)}.kcat;
end
for j=1:n_rxn
    kapp1(j)=rxn{j}.kcat;
end

%% Log Keq's
ln_Keq=log(Keq(jc));
ln_Keq1=[log(Keq1(jc));log(Keq1(jf));log(Keq1(je_v))];

%% Stoichiometric factors
alpha=0.5;
varth(jc)=1;
varth(je_v)=1;
varth(j_fdh1)=alpha;
varth(j_fdh2)=1-alpha;
varth(j_nadh)=2-alpha;
varth(j_nadph)=1+alpha;
varth(j_fdrd)=1;
varth(j_trxrd)=1;
Sc=S(:,jc);
S1=[S(:,jc),S(:,jf),S(:,je_v)];
sc=Sc*varth(1:n_rxn_c);
sf=S(:,jf);
se_v=sum(S(:,je_v),2);
se_a=S(:,je_a);
ss1=sc+sf(:,1)*(2-alpha)+sf(:,2)*(1+alpha)+sf(:,3)+sf(:,4)+se_v+2*se_a;

%% Multi-compartment characteristic concentrations
x_mc_scale=(1-beta_f-beta_e)/(3/x_CO2_Q(i_co2)+2*rxn{j_atps}.nH_transport/x_H_p+5/x_H2_f(i_h2));         % Concentration scale of the multi-compartment reactor system
H2_mc=beta_f*x_c_H2(1)/(beta_f+5*x_mc_scale/x_H2_f(i_h2))/(beta_f+beta_e);                               % H2 concentration in the f stream of the multi-compartment reactor system
x_sc_scale=(1-beta_r)/(3/x_CO2_Q(i_co2)+2*rxn{j_atps}.nH_transport/x_H_p+5/x_H2_f(i_h2));                % Concentration scale of the single-compartment reactor system

%% Reactor optimization subject to crowding constraint without Km (Multi-compartment)
NN=3;
eps_Keq=1e-6;
np_rlx=200;                                             % Number of points used to approximat ...<ln(th0-th1/x) by several ...<a_i*ln(x)+b_i constraints linear in ln(x)
d_rlx=5e-8;                                             % Margin from lower/upper bounds of concentration of oxidized cofactors used to generate points for constraint relaxation
eps_rlx=1e-3;                                           % Margin from the relaxed constraint within which solutions are not accepted. This is to avoid solutions where the hydrogenase ractions of the redox compartment are exactly at equilibrium
n_rxn_cf=n_rxn_c+n_rxn_f*np_rlx;
x_lb=1e-7;
x_lb_atp=1e-3;
x_lb_adp=1e-4;
x_lb_pi=1e-2;
x_lb_ppi=1e-3;
x_ub=0.15;

Aineq=zeros(n_rxn_cf+6,n_met);
bineq=zeros(n_rxn_cf+6,1);
Aeq=zeros(2,n_met);
beq=zeros(2,1);
x_co=zeros(n_rxn_f,np_rlx);
a_rlx=zeros(n_rxn_f,np_rlx);
b_rlx=zeros(n_rxn_f,np_rlx);
theta0=zeros(n_rxn_f,1);
theta1=zeros(n_rxn_f,1);
x_co_lb=zeros(n_rxn_f,1);
lb1=zeros(n_met,1);
ub1=zeros(n_met,1);
alphaE=zeros(n_met,1);
alphaF=zeros(n_met,1);

a1=NN*(NN+1)/2;
a2=NN*(2*NN^2+3*NN+1)/12;

% Constraint relaxation parameters (f compartment)-------------------------
theta0(1)=H2_mc*Keq(jf(1));
theta0(2)=H2_mc*Keq(jf(2));
theta0(3)=sqrt(H2_mc*Keq(jf(3)));
theta0(4)=H2_mc*Keq(jf(4));
theta1(1)=(1+H2_mc*Keq(jf(1)))*(2-alpha)*(beta_f+beta_e)*x_mc_scale/beta_f;
theta1(2)=(1+H2_mc*Keq(jf(2)))*(1+alpha)*(beta_f+beta_e)*x_mc_scale/beta_f;
theta1(3)=(1+sqrt(H2_mc*Keq(jf(3))))*2*x_mc_scale;
theta1(4)=(1+H2_mc*Keq(jf(4)))*(beta_f+beta_e)*x_mc_scale/beta_f;
for i=1:n_rxn_f
    x_co_lb(i)=theta1(i)/theta0(i);
    x_co(i,:)=logspace(log10(x_co_lb(i))+d_rlx,log10(x_ub)-d_rlx,np_rlx);
    for j=1:np_rlx
        a_rlx(i,j)=theta1(i)/(theta0(i)*x_co(i,j)-theta1(i));
        b_rlx(i,j)=log(theta0(i)-theta1(i)/x_co(i,j))-a_rlx(i,j)*log(x_co(i,j));
    end
end
%--------------------------------------------------------------------------

for i=1:n_met
    if i==i_Fdox || i==i_Fdrd
        alphaF(i)=1/beta_f;
        alphaE(i)=0;
    elseif i==i_pyr
        alphaF(i)=(1-eta)/(beta_f+beta_e);
        alphaE(i)=(1-eta)/(beta_f+beta_e);
    else
        alphaF(i)=1/(beta_f+beta_e);
        alphaE(i)=1/(beta_f+beta_e);
    end
end

for i=1:n_met
    if i==i_atp
        lb1(i)=log(x_lb_atp);
    elseif i==i_adp
        lb1(i)=log(x_lb_adp);
    elseif i==i_pi
        lb1(i)=log(x_lb_pi);
    elseif i==i_ppi
        lb1(i)=log(x_lb_ppi);
    elseif i==i_nad
        lb1(i)=log(x_co_lb(1));
    elseif i==i_nadp
        lb1(i)=log(x_co_lb(2));
    elseif i==i_Fdox
        lb1(i)=log(x_co_lb(3));
    elseif i==i_TRXox
        lb1(i)=log(x_co_lb(4));
    elseif i==i_h2
        lb1(i)=-Inf;
    elseif i==i_co2
        lb1(i)=-Inf;
    elseif i==i_Fdrd
        aF=beta_f*alphaF(i)/(beta_f+5*x_mc_scale/x_H2_f(i_h2));
        bF=(sf(i,1)*(2-alpha)+sf(i,2)*(1+alpha)+sf(i,3)+sf(i,4))*x_mc_scale/(beta_f+5*x_mc_scale/x_H2_f(i_h2));
        lb1(i)=log(max([x_lb,x_lb+sc(i)*x_mc_scale,x_lb/alphaF(i),x_lb/aF+bF/aF]));
    else
        aF=beta_f*alphaF(i)/(beta_f+5*x_mc_scale/x_H2_f(i_h2));
        bF=(sf(i,1)*(2-alpha)+sf(i,2)*(1+alpha)+sf(i,3)+sf(i,4))*x_mc_scale/(beta_f+5*x_mc_scale/x_H2_f(i_h2));
        aE=beta_e*alphaE(i)/(beta_e+2*x_mc_scale*rxn{j_atps}.nH_transport/x_H_p);
        bE=(se_v(i)+2*se_a(i))*x_mc_scale/(beta_e+2*x_mc_scale*rxn{j_atps}.nH_transport/x_H_p);
        lb1(i)=log(max([x_lb,x_lb+sc(i)*x_mc_scale,x_lb/alphaF(i),x_lb/alphaE(i),x_lb/aF+bF/aF,x_lb/aE+bE/aE]));
    end
end
for i=1:n_met
    if i==i_h2
        continue;
    elseif i==i_Fdox || i==i_Fdrd
        aF=beta_f*alphaF(i)/(beta_f+5*x_mc_scale/x_H2_f(i_h2));
        bF=(sf(i,1)*(2-alpha)+sf(i,2)*(1+alpha)+sf(i,3)+sf(i,4))*x_mc_scale/(beta_f+5*x_mc_scale/x_H2_f(i_h2));
        ub1(i)=log(min([x_ub,x_ub+sc(i)*x_mc_scale,x_ub/alphaF(i),x_ub/aF+bF/aF]));
    else
        aF=beta_f*alphaF(i)/(beta_f+5*x_mc_scale/x_H2_f(i_h2));
        bF=(sf(i,1)*(2-alpha)+sf(i,2)*(1+alpha)+sf(i,3)+sf(i,4))*x_mc_scale/(beta_f+5*x_mc_scale/x_H2_f(i_h2));
        aE=beta_e*alphaE(i)/(beta_e+2*x_mc_scale*rxn{j_atps}.nH_transport/x_H_p);
        bE=(se_v(i)+2*se_a(i))*x_mc_scale/(beta_e+2*x_mc_scale*rxn{j_atps}.nH_transport/x_H_p);
        ub1(i)=log(min([x_ub,x_ub+sc(i)*x_mc_scale,x_ub/alphaF(i),x_ub/alphaE(i),x_ub/aF+bF/aF,x_ub/aE+bE/aE]));
    end
end

for j=1:n_rxn_c
    w(j)=rxn{jc(j)}.MW*varth(jc(j))/kapp(jc(j));
end
w_=w/norm(w);
W_=diag(w_);
Aineq(1:n_rxn_c,:)=Sc';
ii=1;
for i=1:n_rxn_f
    for j=1:np_rlx
        Aineq(n_rxn_c+ii,i_rd(i))=1;
        Aineq(n_rxn_c+ii,i_ox(i))=-1-a_rlx(i,j);
        bineq(n_rxn_c+ii)=b_rlx(i,j)-eps_rlx;
        ii=ii+1;
    end
end
Aineq(n_rxn_cf+1,i_nadh)=-1;
Aineq(n_rxn_cf+1,i_nad)=1;
Aineq(n_rxn_cf+2,i_nadph)=-1;
Aineq(n_rxn_cf+2,i_nadp)=1;
Aineq(n_rxn_cf+3,i_Fdrd)=-1;
Aineq(n_rxn_cf+3,i_Fdox)=1;
Aineq(n_rxn_cf+4,i_TRXrd)=-1;
Aineq(n_rxn_cf+4,i_TRXox)=1;
Aineq(n_rxn_cf+5,i_atp)=-1;
Aineq(n_rxn_cf+5,i_adp)=1;
Aineq(n_rxn_cf+5,i_pi)=1;
Aineq(n_rxn_cf+6,i_atp)=-1;
Aineq(n_rxn_cf+6,i_amp)=1;
Aineq(n_rxn_cf+6,i_ppi)=1;
Aeq(1,i_pyr)=1;
Aeq(2,i_h2)=1;
bineq(1:n_rxn_c)=ln_Keq-eps_Keq;
beq(1)=log(sc(i_pyr)*x_mc_scale/eta);
beq(2)=log(x_c_H2(1));

Q=a2*Sc*W_*Sc';
c=a1*w_'*Sc'-2*a2*ln_Keq'*W_*Sc';

model.Q=sparse(Q);
model.obj=c;
model.A=sparse([Aineq;Aeq]);
model.rhs=[bineq;beq];
model.lb=lb1;
model.ub=ub1;
model.modelsense='min';
model.sense=strcat(repmat('<',1,n_rxn_cf),'<<<<<<==');
params.NonConvex=2;
sol=gurobi(model,params);
if ~strcmp(sol.status,'OPTIMAL')
    is_mc_feasible=0;
    is_mc_km_feasible=0;
    assimilation_rate=0;
    assimilation_rate_=0;
else 
    is_mc_feasible=1;  
end

if is_mc_feasible
    options=optimoptions('fmincon','Algorithm','active-set','EnableFeasibilityMode',true,'SpecifyObjectiveGradient',true,'SubproblemAlgorithm','cg',...
                    'ScaleProblem',true,'StepTolerance',1e-8,'OptimalityTolerance',5e-8,'ConstraintTolerance',1e-7,'MaxIterations',10000,'MaxFunctionEvaluations',60000);

    [ln_x,~,~,~]=fmincon(@(ln_x)Fun_Obj(ln_x,S,Keq,kapp,varth,rxn,jc,n_rxn_c,n_met),sol.x,Aineq,bineq,Aeq,beq,lb1,ub1,[],options);
    xc=exp(ln_x);

    sum_u=0;
    for j=1:n_rxn_c
        Xi(jc(j))=1;
        for i=1:n_met
            if S(i,jc(j))==0
                continue;
            end
            Xi(jc(j))=Xi(jc(j))*xc(i)^S(i,jc(j));
        end
        Xi(jc(j))=Xi(jc(j))/Keq(jc(j));
        sum_u=sum_u+rxn{jc(j)}.MW*varth(jc(j))/kapp(jc(j))/(1-Xi(jc(j)));
    end
    uc=Enz_tot_v/sum_u;
    
    err=1;
    tol=1e-8;
    damp_f=0.5;
    damp_ev=0.1;
    damp_ea=0.1;
    V_f=1e1;
    V_e=3e1;
    A_e=1e9;
    uf1=(2-alpha)*uc*V_c/V_f;                           
    uf2=(1+alpha)*uc*V_c/V_f;                           
    uf3=uc*V_c/V_f;                                 
    uf4=uc*V_c/V_f;
    ue_v=uc*V_c/V_e; 
    ue_a=2*uc*V_c/A_e;
    fl_Q=-V_c*sc(i_co2)*uc/x_CO2_Q(i_co2);
    fl_H=ue_a*A_e*rxn{j_atps}.nH_transport/x_H_p;
    fl_H2=(uf1+uf2+uf3+uf4)*V_f/x_H2_f(i_h2);
    fl_q=(fl_Q+fl_H+fl_H2)/(1-beta_f-beta_e);
    fl_F=beta_f*fl_q;
    fl_E=beta_e*fl_q;
    while err>tol
        for i=1:n_met
            if i~=i_Fdrd && i~=i_Fdox && i~=i_pyr
                xF(i)=eta_hat*xc(i)/(beta_e+beta_f);
                xE(i)=xF(i);
            elseif i==i_Fdrd || i==i_Fdox && i~=i_pyr
                xF(i)=xc(i)/beta_f;
                xE(i)=0;
            else
                xF(i)=(1-eta)*xc(i)/(beta_e+beta_f);
                xE(i)=xF(i);
            end
            xe(i)=(fl_E*xE(i)+A_e*se_a(i)*ue_a+V_e*se_v(i)*ue_v)/(fl_E+fl_H);
            xf(i)=(fl_F*xF(i)+fl_H2*x_H2_f(i)+V_f*(sf(i,1)*uf1+sf(i,2)*uf2+sf(i,3)*uf3+sf(i,4)*uf4))/(fl_F+fl_H2);
        end
    
        for j=1:n_rxn_f
            Xi(jf(j))=1;
            for i=1:n_met
                if S(i,jf(j))==0
                    continue;
                end
                Xi(jf(j))=Xi(jf(j))*xf(i)^S(i,jf(j));
            end
            Xi(jf(j))=Xi(jf(j))/Keq(jf(j));
        end
        V_f0=V_f;
        V_f=(rxn{jf(1)}.MW*(2-alpha)/kapp(jf(1))/(1-Xi(jf(1)))+rxn{jf(2)}.MW*(1+alpha)/kapp(jf(2))/(1-Xi(jf(2)))+rxn{jf(3)}.MW/kapp(jf(3))/(1-Xi(jf(3)))+rxn{jf(4)}.MW/kapp(jf(4))/(1-Xi(jf(4))))*V_c*uc/Enz_tot_v;
        err_f=abs((V_f-V_f0)/V_f0);
        V_f=damp_f*V_f+(1-damp_f)*V_f0;
    
        for j=1:n_rxn_ev
            Xi(je_v(j))=1;
            for i=1:n_met
                if S(i,je_v(j))==0
                    continue;
                end
                Xi(je_v(j))=Xi(je_v(j))*xe(i)^S(i,je_v(j));
            end
            Xi(je_v(j))=Xi(je_v(j))/Keq(je_v(j));
        end
        V_e0=V_e;
        V_e=(rxn{je_v(1)}.MW/kapp(je_v(1))/(1-Xi(je_v(1)))+rxn{je_v(2)}.MW/kapp(je_v(2))/(1-Xi(je_v(2))))*V_c*uc/Enz_tot_v;
        err_ev=abs((V_e-V_e0)/V_e0);
        V_e=damp_ev*V_e+(1-damp_ev)*V_e0;
    
        Xi(j_atps)=1;
        for i=1:n_met
            if S(i,j_atps)==0
                continue;
            end
            Xi(j_atps)=Xi(j_atps)*xe(i)^S(i,j_atps);
        end
        Xi(j_atps)=Xi(j_atps)/Keq(j_atps);
        kappa_atps=Fun_ATPS_Sat(xe,x_H(1),x_H_p,atps_par,i_atp,i_adp,i_pi);
        A_e0=A_e;
        A_e=2*V_c*uc/kappa_atps/(1-Xi(j_atps))/Enz_tot_a;
        err_ea=abs((A_e-A_e0)/A_e0);
        A_e=damp_ea*A_e+(1-damp_ea)*A_e0;
    
        uf1=(2-alpha)*uc*V_c/V_f;                           
        uf2=(1+alpha)*uc*V_c/V_f;                           
        uf3=uc*V_c/V_f; 
        uf4=uc*V_c/V_f;
        ue_a=2*uc*V_c/A_e;
        fl_Q=-V_c*sc(i_co2)*uc/x_CO2_Q(i_co2);
        fl_H=ue_a*A_e*rxn{j_atps}.nH_transport/x_H_p;
        fl_H2=(uf1+uf2+uf3+uf4)*V_f/x_H2_f(i_h2);
        fl_q=(fl_Q+fl_H+fl_H2)/(1-beta_f-beta_e);
        fl_F=beta_f*fl_q;
        fl_E=beta_e*fl_q;
    
        err=(err_f+err_ev+err_ea)/3;
    end
end

%% Reactor optimization subject to crowding constraint with Km (Multi-compartment)
damp_c=0.2;
w_r=1e-6;
tol=1e-5;
err=1e10;
x0=xc;
xc_(:)=xc(:);
while err>tol && is_mc_feasible
    for j=1:n_rxn_c 
        zeta_(jc(j))=Fun_Zeta(xc_,Km,S,n_met,jc(j));
        w(j)=rxn{jc(j)}.MW*varth(jc(j))/kapp(jc(j))/zeta_(jc(j));
    end
    w_=w/norm(w);
    W_=diag(w_);
    Q=a2*Sc*W_*Sc'+w_r*eye(n_met);
    c=a1*w_'*Sc'-2*a2*ln_Keq'*W_*Sc'-2*w_r*x0';
    model.Q=sparse(Q);
    model.obj=c;
    sol=gurobi(model,params);
    if strcmp(sol.status,'OPTIMAL')
        is_mc_km_feasible=1;
    else 
        is_mc_km_feasible=0;
        break;   
    end
    ln_xc=sol.x;
    xc0=xc_;
    xc_=exp(ln_xc);
    err=norm(xc_-xc0)/n_met;
    xc_=damp_c*xc_+(1-damp_c)*xc0;
end

if is_mc_km_feasible
    options=optimoptions('fmincon','Algorithm','active-set','EnableFeasibilityMode',true,'SpecifyObjectiveGradient',true,'SubproblemAlgorithm','cg',...
                    'ScaleProblem',true,'StepTolerance',1e-8,'OptimalityTolerance',5e-8,'ConstraintTolerance',1e-7,'MaxIterations',10000,'MaxFunctionEvaluations',60000);

    [ln_x,~,~,~]=fmincon(@(ln_x)Fun_Obj_km(ln_x,S,Km,Keq,kapp,varth,rxn,jc,n_rxn_c,n_met),sol.x,Aineq,bineq,Aeq,beq,lb1,ub1,[],options);
    xc_=exp(ln_x);
    
    sum_u=0;
    for j=1:n_rxn_c
        Xi_(jc(j))=1;
        for i=1:n_met
            if S(i,jc(j))==0
                continue;
            end
            Xi_(jc(j))=Xi_(jc(j))*xc_(i)^S(i,jc(j));
        end
        Xi_(jc(j))=Xi_(jc(j))/Keq(jc(j));
        zeta_(jc(j))=Fun_Zeta(xc_,Km,S,n_met,jc(j));
        sum_u=sum_u+rxn{jc(j)}.MW*varth(jc(j))/kapp(jc(j))/(1-Xi_(jc(j)))/zeta_(jc(j));
    end
    uc_=Enz_tot_v/sum_u;

    err=1;
    tol=1e-6;
    damp_f=0.5;
    damp_ev=0.5;
    damp_ea=0.5;
    V_f_=1e0;
    V_e_=1e1;
    A_e_=1e6;
    uf1_=(2-alpha)*uc_*V_c/V_f_;                           
    uf2_=(1+alpha)*uc*V_c/V_f_;                           
    uf3_=uc_*V_c/V_f_;                                 
    uf4_=uc_*V_c/V_f_;
    ue_v_=uc_*V_c/V_e_; 
    ue_a_=2*uc_*V_c/A_e_;
    fl_Q_=-V_c*sc(i_co2)*uc_/x_CO2_Q(i_co2);
    fl_H_=ue_a_*A_e_*rxn{j_atps}.nH_transport/x_H_p;
    fl_H2_=(uf1_+uf2_+uf3_+uf4_)*V_f_/x_H2_f(i_h2);
    fl_q_=(fl_Q_+fl_H_+fl_H2_)/(1-beta_f-beta_e);
    fl_F_=beta_f*fl_q_;
    fl_E_=beta_e*fl_q_;
    while err>tol
        for i=1:n_met
            if i~=i_Fdrd && i~=i_Fdox && i~=i_pyr
                xF_(i)=eta_hat*xc_(i)/(beta_e+beta_f);
                xE_(i)=xF_(i);
            elseif i==i_Fdrd || i==i_Fdox && i~=i_pyr
                xF_(i)=xc_(i)/beta_f;
                xE_(i)=0;
            else
                xF_(i)=(1-eta)*xc_(i)/(beta_e+beta_f);
                xE_(i)=xF_(i);
            end
            xe_(i)=(fl_E_*xE_(i)+A_e_*se_a(i)*ue_a_+V_e_*se_v(i)*ue_v_)/(fl_E_+fl_H_);
            xf_(i)=(fl_F_*xF_(i)+fl_H2_*x_H2_f(i)+V_f_*(sf(i,1)*uf1_+sf(i,2)*uf2_+sf(i,3)*uf3_+sf(i,4)*uf4_))/(fl_F_+fl_H2_);
        end
    
        for j=1:n_rxn_f
            Xi_(jf(j))=1;
            for i=1:n_met
                if S(i,jf(j))==0
                    continue;
                end
                Xi_(jf(j))=Xi_(jf(j))*xf_(i)^S(i,jf(j));
            end
            Xi_(jf(j))=Xi_(jf(j))/Keq(jf(j));
            zeta_(jf(j))=Fun_Zeta(xf_,Km,S,n_met,jf(j));
        end
        V_f0=V_f_;
        V_f_=(rxn{jf(1)}.MW*(2-alpha)/kapp(jf(1))/(1-Xi_(jf(1)))/zeta_(jf(1))+rxn{jf(2)}.MW*(1+alpha)/kapp(jf(2))/(1-Xi_(jf(2)))/zeta_(jf(2))+rxn{jf(3)}.MW/kapp(jf(3))/(1-Xi_(jf(3)))/zeta_(jf(3))+rxn{jf(4)}.MW/kapp(jf(4))/(1-Xi_(jf(4)))/zeta_(jf(4)))*V_c*uc_/Enz_tot_v;
        err_f=abs((V_f_-V_f0)/V_f0);
        V_f_=damp_f*V_f_+(1-damp_f)*V_f0;
    
        for j=1:n_rxn_ev
            Xi_(je_v(j))=1;
            for i=1:n_met
                if S(i,je_v(j))==0
                    continue;
                end
                Xi_(je_v(j))=Xi_(je_v(j))*xe_(i)^S(i,je_v(j));
            end
            Xi_(je_v(j))=Xi_(je_v(j))/Keq(je_v(j));
            zeta_(je_v(j))=Fun_Zeta(xe_,Km,S,n_met,je_v(j));
        end
        V_e0=V_e_;
        V_e_=(rxn{je_v(1)}.MW/kapp(je_v(1))/(1-Xi_(je_v(1)))/zeta_(je_v(1))+rxn{je_v(2)}.MW/kapp(je_v(2))/(1-Xi_(je_v(2)))/zeta_(je_v(2)))*V_c*uc_/Enz_tot_v;
        err_ev=abs((V_e_-V_e0)/V_e0);
        V_e_=damp_ev*V_e_+(1-damp_ev)*V_e0;
    
        Xi_(j_atps)=1;
        for i=1:n_met
            if S(i,j_atps)==0
                continue;
            end
            Xi_(j_atps)=Xi_(j_atps)*xe_(i)^S(i,j_atps);
        end
        Xi_(j_atps)=Xi_(j_atps)/Keq(j_atps);
        kappa_atps_=Fun_ATPS_Sat(xe_,x_H(1),x_H_p,atps_par,i_atp,i_adp,i_pi);
        A_e0=A_e_;
        A_e_=2*V_c*uc_/kappa_atps_/(1-Xi_(j_atps))/Enz_tot_a;
        err_ea=abs((A_e_-A_e0)/A_e0);
        A_e_=damp_ea*A_e_+(1-damp_ea)*A_e0;
    
        uf1_=(2-alpha)*uc_*V_c/V_f_;                           
        uf2_=(1+alpha)*uc_*V_c/V_f_;                           
        uf3_=uc_*V_c/V_f_; 
        uf4_=uc_*V_c/V_f_;
        ue_a_=2*uc_*V_c/A_e_;
        fl_Q_=-V_c*sc(i_co2)*uc_/x_CO2_Q(i_co2);
        fl_H_=ue_a_*A_e_*rxn{j_atps}.nH_transport/x_H_p;
        fl_H2_=(uf1_+uf2_+uf3_+uf4_)*V_f_/x_H2_f(i_h2);
        fl_q_=(fl_Q_+fl_H_+fl_H2_)/(1-beta_f-beta_e);
        fl_F_=beta_f*fl_q_;
        fl_E_=beta_e*fl_q_;
    
        err=(err_f+err_ev+err_ea)/3;
    end
end

%% Reactor optimization subject to crowding constraint without Km (Single-compartment)
eps_Keq=1e-6;
Aineq=zeros(n_rxn1+6,n_met);
bineq=zeros(n_rxn1+6,1);
Aeq=zeros(2,n_met);
beq=zeros(2,1);

for i=1:n_met
    if i==i_atp
        lb1(i)=log(x_lb_atp);
    elseif i==i_adp
        lb1(i)=log(x_lb_adp);
    elseif i==i_pi
        lb1(i)=log(x_lb_pi);
    elseif i==i_ppi
        lb1(i)=log(x_lb_ppi);
    elseif i==i_co2
        lb1(i)=-Inf;
    elseif i==i_h2
        lb1(i)=-Inf;
    else
        lb1(i)=log(x_lb+max(0,ss1(i)*x_sc_scale));
    end
end
for i=1:n_met
    if i~=i_h2
        ub1(i)=log(x_ub+min(0,ss1(i)*x_sc_scale));
    end
end

for j=1:n_rxn_c
    w1(j)=rxn{jc(j)}.MW*varth(jc(j))/kapp1(jc(j));
end
for j=1:n_rxn_f
    w1(n_rxn_c+j)=rxn{jf(j)}.MW*varth(jf(j))/kapp1(jf(j));
end
for j=1:n_rxn_ev
    w1(n_rxn_c+n_rxn_f+j)=rxn{je_v(j)}.MW*varth(je_v(j))/kapp1(je_v(j));
end
w1_=w1/norm(w1);
W1_=diag(w1_);
Aineq(1:n_rxn1,:)=S1';
Aineq(n_rxn1+1,i_nadh)=-1;
Aineq(n_rxn1+1,i_nad)=1;
Aineq(n_rxn1+2,i_nadph)=-1;
Aineq(n_rxn1+2,i_nadp)=1;
Aineq(n_rxn1+3,i_Fdrd)=-1;
Aineq(n_rxn1+3,i_Fdox)=1;
Aineq(n_rxn1+4,i_TRXrd)=-1;
Aineq(n_rxn1+4,i_TRXox)=1;
Aineq(n_rxn1+5,i_atp)=-1;
Aineq(n_rxn1+5,i_adp)=1;
Aineq(n_rxn1+5,i_pi)=1;
Aineq(n_rxn1+6,i_atp)=-1;
Aineq(n_rxn1+6,i_amp)=1;
Aineq(n_rxn1+6,i_ppi)=1;
Aeq(1,i_pyr)=1;
Aeq(2,i_h2)=1;
bineq(1:n_rxn1)=ln_Keq1-eps_Keq;
beq(1)=log(ss1(i_pyr)*x_sc_scale/eta);
beq(2)=log(x_c_H2(2));
Q=a2*S1*W1_*S1';
c=a1*w1_'*S1'-2*a2*ln_Keq1'*W1_*S1';

clear('model');
model.Q=sparse(Q);
model.obj=c;
model.A=sparse([Aineq;Aeq]);
model.rhs=[bineq;beq];
model.modelsense='min';
model.sense=strcat(repmat('<',1,n_rxn1),'<<<<<<==');
model.lb=lb1;
model.ub=ub1;
params.NonConvex=2;
sol=gurobi(model,params);
if ~strcmp(sol.status,'OPTIMAL')
    is_sc_feasible=0;
    is_mc_km_feasible=0;
    assimilation_rate1=0;
    assimilation_rate1_=0;
else 
    is_sc_feasible=1;  

    options=optimoptions('fmincon','Algorithm','active-set','EnableFeasibilityMode',true,'SpecifyObjectiveGradient',true,'SubproblemAlgorithm','cg',...
                    'ScaleProblem',false,'StepTolerance',1e-7,'OptimalityTolerance',1e-5,'ConstraintTolerance',1e-5,'MaxIterations',10000,'MaxFunctionEvaluations',60000);

    [ln_x,~,~,~]=fmincon(@(ln_x)Fun_Obj(ln_x,S,Keq1,kapp1,varth,rxn,[jc;jf;je_v],n_rxn1,n_met),sol.x,Aineq,bineq,Aeq,beq,lb1,ub1,[],options);
    xc1=exp(ln_x);

    sum_u=0;
    for j=1:n_rxn_c
        Xi1(jc(j))=1;
        for i=1:n_met
            if S(i,jc(j))==0
                continue;
            end
            Xi1(jc(j))=Xi1(jc(j))*xc1(i)^S(i,jc(j));
        end
        Xi1(jc(j))=Xi1(jc(j))/Keq1(jc(j));
        sum_u=sum_u+rxn{jc(j)}.MW*varth(jc(j))/kapp1(jc(j))/(1-Xi1(jc(j)));
    end
    for j=1:n_rxn_f
        Xi1(jf(j))=1;
        for i=1:n_met
            if S(i,jf(j))==0
                continue;
            end
            Xi1(jf(j))=Xi1(jf(j))*xc1(i)^S(i,jf(j));
        end
        Xi1(jf(j))=Xi1(jf(j))/Keq1(jf(j));
        sum_u=sum_u+rxn{jf(j)}.MW*varth(jf(j))/kapp1(jf(j))/(1-Xi1(jf(j)));
    end
    for j=1:n_rxn_ev
        Xi1(je_v(j))=1;
        for i=1:n_met
            if S(i,je_v(j))==0
                continue;
            end
            Xi1(je_v(j))=Xi1(je_v(j))*xc1(i)^S(i,je_v(j));
        end
        Xi1(je_v(j))=Xi1(je_v(j))/Keq1(je_v(j));
        sum_u=sum_u+rxn{je_v(j)}.MW*varth(je_v(j))/kapp1(je_v(j))/(1-Xi1(je_v(j)));
    end
    uc1=Enz_tot_v/sum_u;
    
    Xi1(j_atps)=1;
    for i=1:n_met
        if S(i,j_atps)==0
            continue;
        end
        Xi1(j_atps)=Xi1(j_atps)*xc1(i)^S(i,j_atps);
    end
    Xi1(j_atps)=Xi1(j_atps)/Keq1(j_atps);
    kappa_atps1=Fun_ATPS_Sat(xc1,x_H(1),x_H_p,atps_par,i_atp,i_adp,i_pi);
    A_e1=V1*uc1/kappa_atps1/(1-Xi1(j_atps))/Enz_tot_a;
end

%% Reactor optimization subject to crowding constraint with Km (Single-compartment)
damp_c=0.5;
w_r=1e-6;
tol=1e-6;
err=1e10;
x0=xc1;
xc1_(:)=xc1(:);
while err>tol && is_sc_feasible 
    for j=1:n_rxn_c
        zeta1_(jc(j))=Fun_Zeta(xc1_,Km,S,n_met,jc(j));
        w1(j)=rxn{jc(j)}.MW*varth(jc(j))/kapp1(jc(j))/zeta1_(jc(j));
    end
    for j=1:n_rxn_f
        zeta1_(jf(j))=Fun_Zeta(xc1_,Km,S,n_met,jf(j));
        w1(n_rxn_c+j)=rxn{jf(j)}.MW*varth(jf(j))/kapp1(jf(j))/zeta1_(jf(j));
    end
    for j=1:n_rxn_ev
        zeta1_(je_v(j))=Fun_Zeta(xc1_,Km,S,n_met,je_v(j));
        w1(n_rxn_c+n_rxn_f+j)=rxn{je_v(j)}.MW*varth(je_v(j))/kapp1(je_v(j))/zeta1_(je_v(j));
    end
    w1_=w1/norm(w1);
    W1_=diag(w1_);
    Q=a2*S1*W1_*S1'+w_r*eye(n_met);
    c=a1*w1_'*S1'-2*a2*ln_Keq1'*W1_*S1'-2*w_r*x0';
    model.Q=sparse(Q);
    model.obj=c;
    sol=gurobi(model,params);
    if strcmp(sol.status,'OPTIMAL')
        is_sc_km_feasible=1;
    else 
        is_sc_km_feasible=0;
        break;   
    end
    ln_xc=sol.x;
    xc0=xc1_;
    xc1_=exp(ln_xc);
    err=norm(xc1_-xc0)/n_met;
    xc1_=damp_c*xc1_+(1-damp_c)*xc0;
end

if is_sc_km_feasible
    options=optimoptions('fmincon','Algorithm','active-set','EnableFeasibilityMode',true,'SpecifyObjectiveGradient',true,'SubproblemAlgorithm','cg',...
                    'ScaleProblem',false,'StepTolerance',1e-7,'OptimalityTolerance',1e-5,'ConstraintTolerance',1e-5,'MaxIterations',10000,'MaxFunctionEvaluations',60000);

    [ln_x,~,~,~]=fmincon(@(ln_x)Fun_Obj_km(ln_x,S,Km,Keq1,kapp1,varth,rxn,[jc;jf;je_v],n_rxn1,n_met),sol.x,Aineq,bineq,Aeq,beq,lb1,ub1,[],options);
    xc1_=exp(ln_x);

    sum_u=0;
    for j=1:n_rxn_c
        Xi1_(jc(j))=1;
        for i=1:n_met
            if S(i,jc(j))==0
                continue;
            end
            Xi1_(jc(j))=Xi1_(jc(j))*xc1_(i)^S(i,jc(j));
        end
        Xi1_(jc(j))=Xi1_(jc(j))/Keq1(jc(j));
        zeta1_(jc(j))=Fun_Zeta(xc1_,Km,S,n_met,jc(j));
        sum_u=sum_u+rxn{jc(j)}.MW*varth(jc(j))/kapp1(jc(j))/(1-Xi1_(jc(j)))/zeta1_(jc(j));
    end
    for j=1:n_rxn_f
        Xi1_(jf(j))=1;
        for i=1:n_met
            if S(i,jf(j))==0
                continue;
            end
            Xi1_(jf(j))=Xi1_(jf(j))*xc1_(i)^S(i,jf(j));
        end
        Xi1_(jf(j))=Xi1_(jf(j))/Keq1(jf(j));
        zeta1_(jf(j))=Fun_Zeta(xc1_,Km,S,n_met,jf(j));
        sum_u=sum_u+rxn{jf(j)}.MW*varth(jf(j))/kapp1(jf(j))/(1-Xi1_(jf(j)))/zeta1_(jf(j));
    end
    for j=1:n_rxn_ev
        Xi1_(je_v(j))=1;
        for i=1:n_met
            if S(i,je_v(j))==0
                continue;
            end
            Xi1_(je_v(j))=Xi1_(je_v(j))*xc1_(i)^S(i,je_v(j));
        end
        Xi1_(je_v(j))=Xi1_(je_v(j))/Keq1(je_v(j));
        zeta1_(je_v(j))=Fun_Zeta(xc1_,Km,S,n_met,je_v(j));
        sum_u=sum_u+rxn{je_v(j)}.MW*varth(je_v(j))/kapp1(je_v(j))/(1-Xi1_(je_v(j)))/zeta1_(je_v(j));
    end
    uc1_=Enz_tot_v/sum_u;

    Xi1_(j_atps)=1;
    for i=1:n_met
        if S(i,j_atps)==0
            continue;
        end
        Xi1_(j_atps)=Xi1_(j_atps)*xc1_(i)^S(i,j_atps);
    end
    Xi1_(j_atps)=Xi1_(j_atps)/Keq1(j_atps);
    kappa_atps1_=Fun_ATPS_Sat(xc1_,x_H(1),x_H_p,atps_par,i_atp,i_adp,i_pi);
    A_e1_=V1*uc1_/kappa_atps1_/(1-Xi1_(j_atps))/Enz_tot_a;
end

%% Optimal canonical fluxs
if is_mc_feasible
    uf1=(2-alpha)*uc*V_c/V_f;                           % NADH balance
    uf2=(1+alpha)*uc*V_c/V_f;                           % NADPH balance
    uf3=uc*V_c/V_f;                                     % FDrd balance
    uf4=uc*V_c/V_f;                                     % TRXrd balance
    ue_v=uc*V_c/V_e;                                    % AMP, Pi, PPi balance
    ue_a=2*uc*V_c/A_e;                                  % ATP balance
end
if is_mc_km_feasible
    uf1_=(2-alpha)*uc_*V_c/V_f_;                           % NADH balance
    uf2_=(1+alpha)*uc_*V_c/V_f_;                           % NADPH balance
    uf3_=uc_*V_c/V_f_;                                     % FDrd balance
    uf4_=uc_*V_c/V_f_;                                     % TRXrd balance
    ue_v_=uc_*V_c/V_e_;                                    % AMP, Pi, PPi balance
    ue_a_=2*uc_*V_c/A_e_;                                  % ATP balance
end
if is_sc_feasible
    uf1_1=(2-alpha)*uc1;                                % NADH balance               
    uf2_1=(1+alpha)*uc1;                                % NADPH balance                     
    uf3_1=uc1;                                          % FDrd balance                    
    uf4_1=uc1;                                          % TRXrd balance
    ue_v1=uc1;                                          % AMP, Pi, PPi balance
    ue_a1=2*uc1*V1/A_e1;                                % ATP balance
end
if is_sc_km_feasible
    uf1_1_=(2-alpha)*uc1_;                                % NADH balance               
    uf2_1_=(1+alpha)*uc1_;                                % NADPH balance                     
    uf3_1_=uc1_;                                          % FDrd balance                    
    uf4_1_=uc1_;                                          % TRXrd balance
    ue_v1_=uc1_;                                          % AMP, Pi, PPi balance
    ue_a1_=2*uc1_*V1/A_e1_;                                % ATP balance
end

%% Optimal fluxes
if is_mc_feasible
    for j=1:n_rxn_c
        v(jc(j))=varth(jc(j))*uc;
    end
    v(jf(1))=uf1;
    v(jf(2))=uf2;
    v(jf(3))=uf3;
    v(jf(4))=uf4;
    v(je_v(1))=ue_v;
    v(je_v(2))=ue_v;
    v(j_atps)=ue_a;
end
if is_mc_km_feasible
    for j=1:n_rxn_c
        v_(jc(j))=varth(jc(j))*uc_;
    end
    v_(jf(1))=uf1_;
    v_(jf(2))=uf2_;
    v_(jf(3))=uf3_;
    v_(jf(4))=uf4_;
    v_(je_v(1))=ue_v_;
    v_(je_v(2))=ue_v_;
    v_(j_atps)=ue_a_;
end
if is_sc_feasible
    for j=1:n_rxn_c
        v1(jc(j))=varth(jc(j))*uc1;
    end
    v1(jf(1))=uf1_1;
    v1(jf(2))=uf2_1;
    v1(jf(3))=uf3_1;
    v1(jf(4))=uf4_1;
    v1(je_v(1))=ue_v1;
    v1(je_v(2))=ue_v1;
    v1(j_atps)=ue_a1;
end
if is_sc_feasible
    for j=1:n_rxn_c
        v1_(jc(j))=varth(jc(j))*uc1_;
    end
    v1_(jf(1))=uf1_1_;
    v1_(jf(2))=uf2_1_;
    v1_(jf(3))=uf3_1_;
    v1_(jf(4))=uf4_1_;
    v1_(je_v(1))=ue_v1_;
    v1_(je_v(2))=ue_v1_;
    v1_(j_atps)=ue_a1_;
end

%% Optimal flow rates
if is_mc_feasible
    fl_Q=-V_c*sc(i_co2)*uc/x_CO2_Q(i_co2);
    fl_H=ue_a*A_e*rxn{j_atps}.nH_transport/x_H_p;
    fl_H2=(uf1+uf2+uf3+uf4)*V_f/x_H2_f(i_h2);
    fl_q=(fl_Q+fl_H+fl_H2)/(1-beta_f-beta_e);
    fl_F=beta_f*fl_q;
    fl_f=fl_H2+beta_f*fl_q;
    fl_E=beta_e*fl_q;
    fl_e=fl_H+beta_e*fl_q;
    fl_C=(1-beta_f-beta_e)*fl_q;
end
if is_mc_km_feasible
    fl_Q_=-V_c*sc(i_co2)*uc_/x_CO2_Q(i_co2);
    fl_H_=ue_a_*A_e_*rxn{j_atps}.nH_transport/x_H_p;
    fl_H2_=(uf1_+uf2_+uf3_+uf4_)*V_f_/x_H2_f(i_h2);
    fl_q_=(fl_Q_+fl_H_+fl_H2_)/(1-beta_f-beta_e);
    fl_F_=beta_f*fl_q_;
    fl_f_=fl_H2_+beta_f*fl_q_;
    fl_E_=beta_e*fl_q_;
    fl_e_=fl_H_+beta_e*fl_q_;
    fl_C_=(1-beta_f-beta_e)*fl_q_;
end
if is_sc_feasible
    fl_Q_1=-V1*ss1(i_co2)*uc1/x_CO2_Q(i_co2);
    fl_H_1=ue_a1*A_e1*rxn{j_atps}.nH_transport/x_H_p;
    fl_H2_1=(uf1_1+uf2_1+uf3_1+uf4_1)*V1/x_H2_f(i_h2);
    fl_q_1=(fl_Q_1+beta_r*(fl_H_1+fl_H2_1))/(1-beta_r);
    fl_c_1=fl_q_1+fl_H_1+fl_H2_1;
    fl_C_1=(1-beta_r)*fl_c_1;
    fl_R_1=beta_r*fl_c_1;
end
if is_sc_km_feasible
    fl_Q_1_=-V1*ss1(i_co2)*uc1_/x_CO2_Q(i_co2);
    fl_H_1_=ue_a1_*A_e1_*rxn{j_atps}.nH_transport/x_H_p;
    fl_H2_1_=(uf1_1_+uf2_1_+uf3_1_+uf4_1_)*V1/x_H2_f(i_h2);
    fl_q_1_=(fl_Q_1_+beta_r*(fl_H_1_+fl_H2_1_))/(1-beta_r);
    fl_c_1_=fl_q_1_+fl_H_1_+fl_H2_1_;
    fl_C_1_=(1-beta_r)*fl_c_1_;
    fl_R_1_=beta_r*fl_c_1_;
end

%% Optimal metabolite concentrations in q,F,f,E,e,C streams
if is_mc_feasible
    for i=1:n_met
        xq(i)=xc(i)-V_c*sc(i)*uc/fl_q;
    end
    for i=1:n_met
        if i~=i_Fdrd && i~=i_Fdox && i~=i_pyr
            xF(i)=eta_hat*xc(i)/(beta_e+beta_f);
            xE(i)=xF(i);
            xC(i)=(1-eta_hat)*xc(i)/(1-beta_f-beta_e);
        elseif i==i_Fdrd || i==i_Fdox && i~=i_pyr
            xF(i)=xc(i)/beta_f;
            xE(i)=0;
            xC(i)=0;
        else
            xF(i)=(1-eta)*xc(i)/(beta_e+beta_f);
            xE(i)=xF(i);
            xC(i)=eta*xc(i)/(1-beta_f-beta_e);
        end
        xe(i)=(fl_E*xE(i)+A_e*se_a(i)*ue_a+V_e*se_v(i)*ue_v)/(fl_E+fl_H);
        xf(i)=(fl_F*xF(i)+fl_H2*x_H2_f(i)+V_f*(sf(i,1)*uf1+sf(i,2)*uf2+sf(i,3)*uf3+sf(i,4)*uf4))/(fl_F+fl_H2);
    end
end
if is_mc_km_feasible
    for i=1:n_met
        xq_(i)=xc_(i)-V_c*sc(i)*uc_/fl_q_;
    end
    for i=1:n_met
        if i~=i_Fdrd && i~=i_Fdox && i~=i_pyr
            xF_(i)=eta_hat*xc_(i)/(beta_e+beta_f);
            xE_(i)=xF_(i);
            xC_(i)=(1-eta_hat)*xc_(i)/(1-beta_f-beta_e);
        elseif i==i_Fdrd || i==i_Fdox && i~=i_pyr
            xF_(i)=xc_(i)/beta_f;
            xE_(i)=0;
            xC_(i)=0;
        else
            xF_(i)=(1-eta)*xc_(i)/(beta_e+beta_f);
            xE_(i)=xF_(i);
            xC_(i)=eta*xc_(i)/(1-beta_f-beta_e);
        end
        xe_(i)=(fl_E_*xE_(i)+A_e_*se_a(i)*ue_a_+V_e_*se_v(i)*ue_v_)/(fl_E_+fl_H_);
        xf_(i)=(fl_F_*xF_(i)+fl_H2_*x_H2_f(i)+V_f_*(sf(i,1)*uf1_+sf(i,2)*uf2_+sf(i,3)*uf3_+sf(i,4)*uf4_))/(fl_F_+fl_H2_);
    end
end
if is_sc_feasible
    for i=1:n_met
        xq1(i)=(fl_c_1*xc1(i)-fl_H2_1*x_H2_f(i)-V1*ss1(i)*uc1)/fl_q_1;
    end
    for i=1:n_met
        if i==i_pyr
            xR1(i)=(1-eta)*xc1(i)/beta_r;
            xC1(i)=eta*xc1(i)/(1-beta_r);
        else
            xR1(i)=xc1(i)/beta_r;
            xC1(i)=0;
        end
    end
end
if is_sc_km_feasible
    for i=1:n_met
        xq1_(i)=(fl_c_1_*xc1_(i)-fl_H2_1_*x_H2_f(i)-V1*ss1(i)*uc1_)/fl_q_1_;
    end
    for i=1:n_met
        if i==i_pyr
            xR1_(i)=(1-eta)*xc1_(i)/beta_r;
            xC1_(i)=eta*xc1_(i)/(1-beta_r);
        else
            xR1_(i)=xc1_(i)/beta_r;
            xC1_(i)=0;
        end
    end
end

%% Optimal reaction quotients
if is_mc_feasible
    for j=1:n_rxn_c
        Xi(jc(j))=1;
        for i=1:n_met
            if S(i,jc(j))==0
                continue;
            end
            Xi(jc(j))=Xi(jc(j))*xc(i)^S(i,jc(j));
        end
        Xi(jc(j))=Xi(jc(j))/Keq(jc(j));
    end
    for j=1:n_rxn_f
        Xi(jf(j))=1;
        for i=1:n_met
            if S(i,jf(j))==0
                continue;
            end
            Xi(jf(j))=Xi(jf(j))*xf(i)^S(i,jf(j));
        end
        Xi(jf(j))=Xi(jf(j))/Keq(jf(j));
    end
    for j=1:n_rxn_ev
        Xi(je_v(j))=1;
        for i=1:n_met
            if S(i,je_v(j))==0
                continue;
            end
            Xi(je_v(j))=Xi(je_v(j))*xe(i)^S(i,je_v(j));
        end
        Xi(je_v(j))=Xi(je_v(j))/Keq(je_v(j));
    end  
    Xi(j_atps)=1;
    for i=1:n_met
        if S(i,j_atps)==0
            continue;
        end
        Xi(j_atps)=Xi(j_atps)*xe(i)^S(i,j_atps);
    end
    Xi(j_atps)=Xi(j_atps)/Keq(j_atps);
end
if is_mc_km_feasible
    for j=1:n_rxn_c
        Xi_(jc(j))=1;
        for i=1:n_met
            if S(i,jc(j))==0
                continue;
            end
            Xi_(jc(j))=Xi_(jc(j))*xc_(i)^S(i,jc(j));
        end
        Xi_(jc(j))=Xi_(jc(j))/Keq(jc(j));
    end
    for j=1:n_rxn_f
        Xi_(jf(j))=1;
        for i=1:n_met
            if S(i,jf(j))==0
                continue;
            end
            Xi_(jf(j))=Xi_(jf(j))*xf_(i)^S(i,jf(j));
        end
        Xi_(jf(j))=Xi_(jf(j))/Keq(jf(j));
    end
    for j=1:n_rxn_ev
        Xi_(je_v(j))=1;
        for i=1:n_met
            if S(i,je_v(j))==0
                continue;
            end
            Xi_(je_v(j))=Xi_(je_v(j))*xe_(i)^S(i,je_v(j));
        end
        Xi_(je_v(j))=Xi_(je_v(j))/Keq(je_v(j));
    end  
    Xi_(j_atps)=1;
    for i=1:n_met
        if S(i,j_atps)==0
            continue;
        end
        Xi_(j_atps)=Xi_(j_atps)*xe_(i)^S(i,j_atps);
    end
    Xi_(j_atps)=Xi_(j_atps)/Keq(j_atps);
end
if is_sc_feasible
    for j=1:n_rxn
        Xi1(j)=1;
        for i=1:n_met
            if S(i,j)==0
                continue;
            end
            Xi1(j)=Xi1(j)*xc1(i)^S(i,j);
        end
        Xi1(j)=Xi1(j)/Keq1(j);
    end
end
if is_sc_km_feasible
    for j=1:n_rxn
        Xi1_(j)=1;
        for i=1:n_met
            if S(i,j)==0
                continue;
            end
            Xi1_(j)=Xi1_(j)*xc1_(i)^S(i,j);
        end
        Xi1_(j)=Xi1_(j)/Keq1(j);
    end
end

%% Optimal enzyme concentrations
if is_mc_feasible
    for j=1:n_rxn_c
        Enz(jc(j))=varth(jc(j))*uc/kapp(jc(j))/(1-Xi(jc(j)));
    end
    Enz(jf(1))=uf1/kapp(jf(1))/(1-Xi(jf(1)));
    Enz(jf(2))=uf2/kapp(jf(2))/(1-Xi(jf(2)));
    Enz(jf(3))=uf3/kapp(jf(3))/(1-Xi(jf(3)));
    Enz(jf(4))=uf4/kapp(jf(4))/(1-Xi(jf(4)));
    Enz(je_v(1))=ue_v/kapp(je_v(1))/(1-Xi(je_v(1)));
    Enz(je_v(2))=ue_v/kapp(je_v(2))/(1-Xi(je_v(2)));
    Enz(j_atps)=ue_a/kappa_atps/(1-Xi(j_atps));
end
if is_mc_km_feasible
    for j=1:n_rxn_c
        Enz_(jc(j))=varth(jc(j))*uc_/kapp(jc(j))/(1-Xi_(jc(j)))/zeta_(jc(j));
    end
    Enz_(jf(1))=uf1_/kapp(jf(1))/(1-Xi_(jf(1)))/zeta_(jf(1));
    Enz_(jf(2))=uf2_/kapp(jf(2))/(1-Xi_(jf(2)))/zeta_(jf(2));
    Enz_(jf(3))=uf3_/kapp(jf(3))/(1-Xi_(jf(3)))/zeta_(jf(3));
    Enz_(jf(4))=uf4_/kapp(jf(4))/(1-Xi_(jf(4)))/zeta_(jf(4));
    Enz_(je_v(1))=ue_v_/kapp(je_v(1))/(1-Xi_(je_v(1)))/zeta_(je_v(1));
    Enz_(je_v(2))=ue_v_/kapp(je_v(2))/(1-Xi_(je_v(2)))/zeta_(je_v(2));
    Enz_(j_atps)=ue_a_/kappa_atps_/(1-Xi_(j_atps));
end
if is_sc_feasible
    for j=1:n_rxn_c
        Enz1(jc(j))=varth(jc(j))*uc1/kapp1(jc(j))/(1-Xi1(jc(j)));
    end
    Enz1(jf(1))=uf1_1/kapp1(jf(1))/(1-Xi1(jf(1)));
    Enz1(jf(2))=uf2_1/kapp1(jf(2))/(1-Xi1(jf(2)));
    Enz1(jf(3))=uf3_1/kapp1(jf(3))/(1-Xi1(jf(3)));
    Enz1(jf(4))=uf4_1/kapp1(jf(4))/(1-Xi1(jf(4)));
    Enz1(je_v(1))=ue_v1/kapp1(je_v(1))/(1-Xi1(je_v(1)));
    Enz1(je_v(2))=ue_v1/kapp1(je_v(2))/(1-Xi1(je_v(2)));
    Enz1(j_atps)=ue_a1/kappa_atps1/(1-Xi1(j_atps));
end
if is_sc_km_feasible
    for j=1:n_rxn_c
        Enz1_(jc(j))=varth(jc(j))*uc1_/kapp1(jc(j))/(1-Xi1_(jc(j)))/zeta1_(jc(j));
    end
    Enz1_(jf(1))=uf1_1_/kapp1(jf(1))/(1-Xi1_(jf(1)))/zeta1_(jf(1));
    Enz1_(jf(2))=uf2_1_/kapp1(jf(2))/(1-Xi1_(jf(2)))/zeta1_(jf(2));
    Enz1_(jf(3))=uf3_1_/kapp1(jf(3))/(1-Xi1_(jf(3)))/zeta1_(jf(3));
    Enz1_(jf(4))=uf4_1_/kapp1(jf(4))/(1-Xi1_(jf(4)))/zeta1_(jf(4));
    Enz1_(je_v(1))=ue_v1_/kapp1(je_v(1))/(1-Xi1_(je_v(1)))/zeta1_(je_v(1));
    Enz1_(je_v(2))=ue_v1_/kapp1(je_v(2))/(1-Xi1_(je_v(2)))/zeta1_(je_v(2));
    Enz1_(j_atps)=ue_a1_/kappa_atps1_/(1-Xi1_(j_atps));
end

%% Optimal partition of Carbon compartment
if is_mc_feasible
    sum_f=0;
    for j=1:size(jc_f,1)
        sum_f=sum_f+rxn{jc_f(j)}.MW*Enz(jc_f(j));
    end
    sum_m=0;
    for j=1:size(jc_m,1)
        sum_m=sum_m+rxn{jc_m(j)}.MW*Enz(jc_m(j));
    end
    r_fm=sum_f/sum_m;
    phi_opt=r_fm/(1+r_fm);
end
if is_mc_km_feasible
    sum_f=0;
    for j=1:size(jc_f,1)
        sum_f=sum_f+rxn{jc_f(j)}.MW*Enz_(jc_f(j));
    end
    sum_m=0;
    for j=1:size(jc_m,1)
        sum_m=sum_m+rxn{jc_m(j)}.MW*Enz_(jc_m(j));
    end
    r_fm=sum_f/sum_m;
    phi_opt_=r_fm/(1+r_fm);
end

%% Assimilation rate
if is_mc_feasible
    assimilation_rate=fl_Q*x_CO2_Q(i_co2);
end
if is_mc_km_feasible
    assimilation_rate_=fl_Q_*x_CO2_Q(i_co2);
end
if is_sc_feasible
    assimilation_rate1=fl_Q_1*x_CO2_Q(i_co2);
end
if is_sc_km_feasible
    assimilation_rate1_=fl_Q_1_*x_CO2_Q(i_co2);
end

%% Checking mass balances
tol_mb=1e-8;
if is_mc_feasible
    err_mc_in=fl_Q*x_CO2_Q+fl_f*xf+fl_e*xe-fl_q*xq;
    err_mc_out=fl_q*xc-fl_F*xF-fl_F*xE-fl_C*xC;
    is_mc_in_mb=all(abs(err_mc_in)<tol_mb);
    is_mc_out_mb=all(abs(err_mc_out)<tol_mb);
end
if is_sc_feasible
    err_sc_in=fl_Q_1*x_CO2_Q+fl_R_1*xR1-fl_q_1*xq1;
    err_sc_out=fl_c_1*xc1-fl_R_1*xR1-fl_C_1*xC1;
    is_sc_in_mb=all(abs(err_sc_in)<tol_mb);
    is_sc_out_mb=all(abs(err_sc_out)<tol_mb);
end

%% Store solutions
if is_mc_feasible
    multi_comp.feas=1;
    multi_comp.xq=xq;
    multi_comp.xF=xF;
    multi_comp.xf=xf;
    multi_comp.xE=xE;
    multi_comp.xe=xe;
    multi_comp.xC=xC;
    multi_comp.xc=xc;
    multi_comp.Enz=Enz;
    multi_comp.v=v;
    multi_comp.fl_Q=fl_Q;
    multi_comp.fl_H=fl_H;
    multi_comp.fl_H2=fl_H2;
    multi_comp.fl_q=fl_q;
    multi_comp.fl_F=fl_F;
    multi_comp.fl_f=fl_f;
    multi_comp.fl_E=fl_E;
    multi_comp.fl_e=fl_e;
    multi_comp.fl_C=fl_C;
    multi_comp.assimilation_rate=assimilation_rate;
    multi_comp.phi_opt=phi_opt;
    multi_comp.kappa_atps=kappa_atps;
else
    multi_comp.feas=0;
    multi_comp.xq=zeros(n_met,1);
    multi_comp.xF=zeros(n_met,1);
    multi_comp.xf=zeros(n_met,1);
    multi_comp.xE=zeros(n_met,1);
    multi_comp.xe=zeros(n_met,1);
    multi_comp.xC=zeros(n_met,1);
    multi_comp.xc=zeros(n_met,1);
    multi_comp.Enz=zeros(n_rxn,1);
    multi_comp.v=zeros(n_rxn,1);
    multi_comp.fl_Q=0;
    multi_comp.fl_H=0;
    multi_comp.fl_H2=0;
    multi_comp.fl_q=0;
    multi_comp.fl_F=0;
    multi_comp.fl_f=0;
    multi_comp.fl_E=0;
    multi_comp.fl_e=0;
    multi_comp.fl_C=0;
    multi_comp.assimilation_rate=0;
    multi_comp.phi_opt=0;
    multi_comp.kappa_atps=0;
end
if is_mc_km_feasible
    multi_comp.feas_=1;
    multi_comp.xq_=xq_;
    multi_comp.xF_=xF_;
    multi_comp.xf_=xf_;
    multi_comp.xE_=xE_;
    multi_comp.xe_=xe_;
    multi_comp.xC_=xC_;
    multi_comp.xc_=xc_;
    multi_comp.Enz_=Enz_;
    multi_comp.v_=v_;
    multi_comp.fl_Q_=fl_Q_;
    multi_comp.fl_H_=fl_H_;
    multi_comp.fl_H2_=fl_H2_;
    multi_comp.fl_q_=fl_q_;
    multi_comp.fl_F_=fl_F_;
    multi_comp.fl_f_=fl_f_;
    multi_comp.fl_E_=fl_E_;
    multi_comp.fl_e_=fl_e_;
    multi_comp.fl_C_=fl_C_;
    multi_comp.assimilation_rate_=assimilation_rate_;
    multi_comp.phi_opt_=phi_opt_;
    multi_comp.kappa_atps_=kappa_atps_;
else
    multi_comp.feas_=0;
    multi_comp.xq_=zeros(n_met,1);
    multi_comp.xF_=zeros(n_met,1);
    multi_comp.xf_=zeros(n_met,1);
    multi_comp.xE_=zeros(n_met,1);
    multi_comp.xe_=zeros(n_met,1);
    multi_comp.xC_=zeros(n_met,1);
    multi_comp.xc_=zeros(n_met,1);
    multi_comp.Enz_=zeros(n_rxn,1);
    multi_comp.v_=zeros(n_rxn,1);
    multi_comp.fl_Q_=0;
    multi_comp.fl_H_=0;
    multi_comp.fl_H2_=0;
    multi_comp.fl_q_=0;
    multi_comp.fl_F_=0;
    multi_comp.fl_f_=0;
    multi_comp.fl_E_=0;
    multi_comp.fl_e_=0;
    multi_comp.fl_C_=0;
    multi_comp.assimilation_rate_=0;
    multi_comp.phi_opt_=0;
    multi_comp.kappa_atps_=0;
end
if is_sc_feasible
    single_comp.feas=1;
    single_comp.xq=xq1;
    single_comp.xC=xC1;
    single_comp.xc=xc1;
    single_comp.xR=xR1;
    single_comp.Enz=Enz1;
    single_comp.v=v1;
    single_comp.fl_Q=fl_Q_1;
    single_comp.fl_H=fl_H_1;
    single_comp.fl_H2=fl_H2_1;
    single_comp.fl_q=fl_q_1;
    single_comp.fl_c=fl_c_1;
    single_comp.fl_C=fl_C_1;
    single_comp.fl_R=fl_R_1;
    single_comp.assimilation_rate=assimilation_rate1;
    single_comp.kappa_atps=kappa_atps1;
else
    single_comp.feas=0;
    single_comp.xq=zeros(n_met,1);
    single_comp.xC=zeros(n_met,1);
    single_comp.xc=zeros(n_met,1);
    single_comp.xR=zeros(n_met,1);
    single_comp.Enz=zeros(n_rxn,1);
    single_comp.v=zeros(n_rxn,1);
    single_comp.fl_Q=0;
    single_comp.fl_H=0;
    single_comp.fl_H2=0;
    single_comp.fl_q=0;
    single_comp.fl_c=0;
    single_comp.fl_C=0;
    single_comp.fl_R=0;
    single_comp.assimilation_rate=0;
    single_comp.kappa_atps=0;
end
if is_sc_km_feasible
    single_comp.feas_=1;
    single_comp.xq_=xq1_;
    single_comp.xC_=xC1_;
    single_comp.xc_=xc1_;
    single_comp.xR_=xR1_;
    single_comp.Enz_=Enz1_;
    single_comp.v_=v1_;
    single_comp.fl_Q_=fl_Q_1_;
    single_comp.fl_H_=fl_H_1_;
    single_comp.fl_H2_=fl_H2_1_;
    single_comp.fl_q_=fl_q_1_;
    single_comp.fl_c_=fl_c_1_;
    single_comp.fl_C_=fl_C_1_;
    single_comp.fl_R_=fl_R_1_;
    single_comp.assimilation_rate_=assimilation_rate1_;
    single_comp.kappa_atps_=kappa_atps1_;
else
    single_comp.feas_=0;
    single_comp.xq_=zeros(n_met,1);
    single_comp.xC_=zeros(n_met,1);
    single_comp.xc_=zeros(n_met,1);
    single_comp.xR_=zeros(n_met,1);
    single_comp.Enz_=zeros(n_rxn,1);
    single_comp.v_=zeros(n_rxn,1);
    single_comp.fl_Q_=0;
    single_comp.fl_H_=0;
    single_comp.fl_H2_=0;
    single_comp.fl_q_=0;
    single_comp.fl_c_=0;
    single_comp.fl_C_=0;
    single_comp.fl_R_=0;
    single_comp.assimilation_rate_=0;
    single_comp.kappa_atps_=0;
end

%% Save MAT
save(strSaveFileMAT,'multi_comp','single_comp','Xi','Xi1','Xi_','Xi1_','zeta_','zeta1_','Keq','Keq1','kapp','kapp1');

%% Display
if is_mc_feasible
    fprintf('Multi-compartment solution: ------------------------------- \n');
    fprintf('Assimilation Rate (mol/s) = %e\n',assimilation_rate);
    fprintf('Carbon Module Flux (M/s) = %e\n',uc);
    fprintf('Pyruvate Outlet Concentration (M) = %e\n',xC(i_pyr));
    fprintf('Volume of Redox Module (L) = %e\n',V_f);
    fprintf('Area of Energy Module (m^2) = %e\n',A_e);
    fprintf('Optimal phi = %e\n',phi_opt);
    fprintf('fl_Q (L/s) = %e\n',fl_Q);
    fprintf('fl_q (L/s) = %e\n',fl_q);
    fprintf('fl_F (L/s) = %e\n',fl_F);
    fprintf('fl_f (L/s) = %e\n',fl_f);
    fprintf('fl_E (L/s) = %e\n',fl_E);
    fprintf('fl_e (L/s) = %e\n',fl_e);
    fprintf('fl_C (L/s) = %e\n',fl_C);
    fprintf('fl_c (L/s) = %e\n',fl_q);
    fprintf('fl_H (L/s) = %e\n',fl_H);
    fprintf('fl_H2 (L/s) = %e\n',fl_H2);
    fprintf('\n');
    if is_mc_in_mb
        fprintf('Mass balance is satisfied at the inlet node!\n');
    else
        fprintf('Mass balance is not quite satisfied at the inlet node .... hmmmm\n');
    end
    if is_mc_out_mb
        fprintf('Mass balance is satisfied at the outlet node!\n');
    else
        fprintf('Mass balance is not quite satisfied at the outlet node .... hmmmm\n');
    end
else
    fprintf('Infeasiblily encountered when optimizing the multi-compartment configuration.\n');
end
if is_sc_feasible
    fprintf('\n');
    fprintf('Single-compartment solution: ------------------------------- \n');
    fprintf('Assimilation Rate (mol/s) = %e\n',assimilation_rate1);
    fprintf('Carbon Module Flux (M/s) = %e\n',uc1);
    fprintf('Pyruvate Outlet Concentration (M) = %e\n',xC1(i_pyr));
    fprintf('Area of Energy Module (m^2) = %e\n',A_e1);
    fprintf('fl_Q (L/s) = %e\n',fl_Q_1);
    fprintf('fl_q (L/s) = %e\n',fl_q_1);
    fprintf('fl_C (L/s) = %e\n',fl_C_1);
    fprintf('fl_c (L/s) = %e\n',fl_c_1);
    fprintf('fl_R (L/s) = %e\n',fl_R_1);
    fprintf('fl_H (L/s) = %e\n',fl_H_1);
    fprintf('fl_H2 (L/s) = %e\n',fl_H2_1);
    fprintf('\n');
    if is_sc_in_mb
        fprintf('Mass balance is satisfied at the inlet node!\n');
    else
        fprintf('Mass balance is not quite satisfied at the inlet node .... hmmmm\n');
    end
    if is_sc_out_mb
        fprintf('Mass balance is satisfied at the outlet node!\n');
    else
        fprintf('Mass balance is not quite satisfied at the outlet node .... hmmmm\n');
    end
else
    fprintf('Infeasiblily encountered when optimizing the single-compartment configuration.\n');
end

%% Plot specifications 
Enz_BaseValue=1e-15;
Enz_right=Enz_BaseValue*ones(n_rxn,1);
Enz_right(j_atps)=Enz(j_atps);
Enz_right1=Enz_BaseValue*ones(n_rxn,1);
Enz_right1(j_atps)=Enz1(j_atps);
Enz_right_=Enz_BaseValue*ones(n_rxn,1);
Enz_right_(j_atps)=Enz_(j_atps);
Enz_right1_=Enz_BaseValue*ones(n_rxn,1);
Enz_right1_(j_atps)=Enz1_(j_atps);
Enz_color_left=[156, 23, 94]/255;
Enz_color_right=[30, 109, 212]/255;
Enz_color_left1=[255, 161, 211]/255;
Enz_color_right1=[143, 192, 255]/255;
Enz_color_left_=[68, 50, 168]/255;
Enz_color_right_=[57, 163, 147]/255;
Enz_color_left1_=[147, 130, 237]/255;
Enz_color_right1_=[132, 186, 178]/255;

v_BaseValue=1e-9;
v_right=v_BaseValue*ones(n_rxn,1);
v_right(j_atps)=v(j_atps);
v_right1=v_BaseValue*ones(n_rxn,1);
v_right1(j_atps)=v1(j_atps);
v_right_=v_BaseValue*ones(n_rxn,1);
v_right_(j_atps)=v_(j_atps);
v_right1_=v_BaseValue*ones(n_rxn,1);
v_right1_(j_atps)=v1_(j_atps);
v_color_left=[156, 23, 94]/255;
v_color_right=[30, 109, 212]/255;
v_color_left1=[255, 161, 211]/255;
v_color_right1=[143, 192, 255]/255;
v_color_left_=[68, 50, 168]/255;
v_color_right_=[57, 163, 147]/255;
v_color_left1_=[147, 130, 237]/255;
v_color_right1_=[132, 186, 178]/255;

x_BaseValue=0.1*x_lb;
xq_color=[255, 164, 31]/255;
xc_color=[50, 178, 50]/255;

BaseValue_Ar=1e-4;
assim_rate=[assimilation_rate_ assimilation_rate-assimilation_rate_;assimilation_rate1_ assimilation_rate1-assimilation_rate1_];
ar_color_bottom=[33, 33, 33]/255;
ar_color_top=[255, 200, 3]/255;
x_label={'MC';'SC'};

%% Plots
if is_mc_feasible==1 && is_sc_feasible==1
    hh1=figure;
    mu_factor=1e6;
    yyaxis left
    bar1=bar([Enz,Enz1]*mu_factor);
    ylim([mu_factor*Enz_BaseValue,1e6]);
    set(gca,'XTick',1:n_rxn,'XTickLabel',rxn_id);
    xtickangle(gca,90);
    ylabel('$E$ ($\mu$M)','FontSize',20,'FontName','Arial','Interpreter','latex');
    set(bar1(1),'BaseValue',mu_factor*Enz_BaseValue,'FaceColor',Enz_color_left,'EdgeColor','black');
    set(bar1(2),'BaseValue',mu_factor*Enz_BaseValue,'FaceColor',Enz_color_left1,'EdgeColor','black');
    set(gca,'YScale','log','YColor',Enz_color_left,'FontSize',20,'FontName','Arial');
    set(gcf,'PaperUnits','centimeters');
    set(gcf,'PaperSize',[20 12]);
    set(gcf,'PaperPosition',[0 0 20 12]);
    
    yyaxis right
    bar2=bar([Enz_right,Enz_right1]*mu_factor);
    ylim([mu_factor*Enz_BaseValue,1e-2]);
    ylabel('$E$ ($\mu$mol/m$^2$)','FontSize',20,'FontName','Arial','Interpreter','latex');
    set(bar2(1),'BaseValue',mu_factor*Enz_BaseValue,'FaceColor',Enz_color_right,'EdgeColor','black');
    set(bar2(2),'BaseValue',mu_factor*Enz_BaseValue,'FaceColor',Enz_color_right1,'EdgeColor','black');
    set(gca,'YScale','log','YColor',Enz_color_right,'FontSize',20,'FontName','Arial');
    
    hh2=figure;
    mu_factor=1e0;
    yyaxis left
    bar1=bar([v,v1]*mu_factor);
    ylim([mu_factor*v_BaseValue,1e1]);
    set(gca,'XTick',1:n_rxn,'XTickLabel',rxn_id);
    xtickangle(gca,90);
    ylabel('$v$ (M/s)','FontSize',20,'FontName','Arial','Interpreter','latex');
    set(bar1(1),'BaseValue',mu_factor*v_BaseValue,'FaceColor',v_color_left,'EdgeColor','black');
    set(bar1(2),'BaseValue',mu_factor*v_BaseValue,'FaceColor',v_color_left1,'EdgeColor','black');
    set(gca,'YScale','log','YColor',v_color_left,'FontSize',20,'FontName','Arial');
    set(gcf,'PaperUnits','centimeters');
    set(gcf,'PaperSize',[20 12]);
    set(gcf,'PaperPosition',[0 0 20 12]);
    
    yyaxis right
    bar2=bar([v_right,v_right1]*mu_factor);
    ylim([mu_factor*v_BaseValue,1e-7]);
    ylabel('$v$ (mol/m$^2/s$)','FontSize',20,'FontName','Arial','Interpreter','latex');
    set(bar2(1),'BaseValue',mu_factor*Enz_BaseValue,'FaceColor',v_color_right,'EdgeColor','black');
    set(bar2(2),'BaseValue',mu_factor*Enz_BaseValue,'FaceColor',v_color_right1,'EdgeColor','black');
    set(gca,'YScale','log','YColor',v_color_right,'FontSize',20,'FontName','Arial');
    
    hh3=figure;
    bar1=bar([xq,xc]);
    ylim([x_BaseValue,1e0]);
    set(gca,'XTick',1:n_met,'XTickLabel',met_id);
    legend('q-stream','c-stream');
    xtickangle(gca,90);
    ylabel('$x$ (M)','FontSize',20,'FontName','Arial','Interpreter','latex');
    set(bar1(1),'BaseValue',x_BaseValue,'FaceColor',xq_color,'EdgeColor','black');
    set(bar1(2),'BaseValue',x_BaseValue,'FaceColor',xc_color,'EdgeColor','black');
    set(gca,'YScale','log','FontSize',20,'FontName','Arial');
    set(gcf,'PaperUnits','centimeters');
    set(gcf,'PaperSize',[25 12]);
    set(gcf,'PaperPosition',[0 0 25 12]);
end
if is_mc_km_feasible==1 && is_sc_km_feasible==1
    hh4=figure;
    mu_factor=1e6;
    yyaxis left
    bar1=bar([Enz_,Enz1_]*mu_factor);
    ylim([mu_factor*Enz_BaseValue,1e6]);
    set(gca,'XTick',1:n_rxn,'XTickLabel',rxn_id);
    xtickangle(gca,90);
    ylabel('$E$ ($\mu$M)','FontSize',20,'FontName','Arial','Interpreter','latex');
    set(bar1(1),'BaseValue',mu_factor*Enz_BaseValue,'FaceColor',Enz_color_left_,'EdgeColor','black');
    set(bar1(2),'BaseValue',mu_factor*Enz_BaseValue,'FaceColor',Enz_color_left1_,'EdgeColor','black');
    set(gca,'YScale','log','YColor',Enz_color_left_,'FontSize',20,'FontName','Arial');
    set(gcf,'PaperUnits','centimeters');
    set(gcf,'PaperSize',[20 12]);
    set(gcf,'PaperPosition',[0 0 20 12]);
    
    yyaxis right
    bar2=bar([Enz_right_,Enz_right1_]*mu_factor);
    ylim([mu_factor*Enz_BaseValue,1e-2]);
    ylabel('$E$ ($\mu$mol/m$^2$)','FontSize',20,'FontName','Arial','Interpreter','latex');
    set(bar2(1),'BaseValue',mu_factor*Enz_BaseValue,'FaceColor',Enz_color_right_,'EdgeColor','black');
    set(bar2(2),'BaseValue',mu_factor*Enz_BaseValue,'FaceColor',Enz_color_right1_,'EdgeColor','black');
    set(gca,'YScale','log','YColor',Enz_color_right_,'FontSize',20,'FontName','Arial');
    
    hh5=figure;
    mu_factor=1e0;
    yyaxis left
    bar1=bar([v_,v1_]*mu_factor);
    ylim([mu_factor*v_BaseValue,1e1]);
    set(gca,'XTick',1:n_rxn,'XTickLabel',rxn_id);
    xtickangle(gca,90);
    ylabel('$v$ (M/s)','FontSize',20,'FontName','Arial','Interpreter','latex');
    set(bar1(1),'BaseValue',mu_factor*v_BaseValue,'FaceColor',v_color_left_,'EdgeColor','black');
    set(bar1(2),'BaseValue',mu_factor*v_BaseValue,'FaceColor',v_color_left1_,'EdgeColor','black');
    set(gca,'YScale','log','YColor',v_color_left_,'FontSize',20,'FontName','Arial');
    set(gcf,'PaperUnits','centimeters');
    set(gcf,'PaperSize',[20 12]);
    set(gcf,'PaperPosition',[0 0 20 12]);
    
    yyaxis right
    bar2=bar([v_right_,v_right1_]*mu_factor);
    ylim([mu_factor*v_BaseValue,1e-7]);
    ylabel('$v$ (mol/m$^2/s$)','FontSize',20,'FontName','Arial','Interpreter','latex');
    set(bar2(1),'BaseValue',mu_factor*Enz_BaseValue,'FaceColor',v_color_right_,'EdgeColor','black');
    set(bar2(2),'BaseValue',mu_factor*Enz_BaseValue,'FaceColor',v_color_right1_,'EdgeColor','black');
    set(gca,'YScale','log','YColor',v_color_right_,'FontSize',20,'FontName','Arial');
    
    hh6=figure;
    bar1=bar([xq_,xc_]);
    ylim([x_BaseValue,1e0]);
    set(gca,'XTick',1:n_met,'XTickLabel',met_id);
    legend('q-stream','c-stream');
    xtickangle(gca,90);
    ylabel('$x$ (M)','FontSize',20,'FontName','Arial','Interpreter','latex');
    set(bar1(1),'BaseValue',x_BaseValue,'FaceColor',xq_color,'EdgeColor','black');
    set(bar1(2),'BaseValue',x_BaseValue,'FaceColor',xc_color,'EdgeColor','black');
    set(gca,'YScale','log','FontSize',20,'FontName','Arial');
    set(gcf,'PaperUnits','centimeters');
    set(gcf,'PaperSize',[25 12]);
    set(gcf,'PaperPosition',[0 0 25 12]);
end
if is_mc_feasible==1
    hh7=figure;
    mu_factor=1e6;
    yyaxis left
    bar1=bar([Enz,Enz_]*mu_factor);
    ylim([mu_factor*Enz_BaseValue,1e6]);
    set(gca,'XTick',1:n_rxn,'XTickLabel',rxn_id);
    xtickangle(gca,90);
    ylabel('$E$ ($\mu$M)','FontSize',20,'FontName','Arial','Interpreter','latex');
    set(bar1(1),'BaseValue',mu_factor*Enz_BaseValue,'FaceColor',Enz_color_left,'EdgeColor','black');
    set(bar1(2),'BaseValue',mu_factor*Enz_BaseValue,'FaceColor',Enz_color_left_,'EdgeColor','black');
    set(gca,'YScale','log','YColor',Enz_color_left,'FontSize',20,'FontName','Arial');
    set(gcf,'PaperUnits','centimeters');
    set(gcf,'PaperSize',[20 12]);
    set(gcf,'PaperPosition',[0 0 20 12]);
    
    yyaxis right
    bar2=bar([Enz_right,Enz_right_]*mu_factor);
    ylim([mu_factor*Enz_BaseValue,1e-2]);
    ylabel('$E$ ($\mu$mol/m$^2$)','FontSize',20,'FontName','Arial','Interpreter','latex');
    set(bar2(1),'BaseValue',mu_factor*Enz_BaseValue,'FaceColor',Enz_color_right,'EdgeColor','black');
    set(bar2(2),'BaseValue',mu_factor*Enz_BaseValue,'FaceColor',Enz_color_right_,'EdgeColor','black');
    set(gca,'YScale','log','YColor',Enz_color_right,'FontSize',20,'FontName','Arial');
end
if is_sc_feasible==1 
    hh8=figure;
    mu_factor=1e6;
    yyaxis left
    bar1=bar([Enz1,Enz1_]*mu_factor);
    ylim([mu_factor*Enz_BaseValue,1e6]);
    set(gca,'XTick',1:n_rxn,'XTickLabel',rxn_id);
    xtickangle(gca,90);
    ylabel('$E$ ($\mu$M)','FontSize',20,'FontName','Arial','Interpreter','latex');
    set(bar1(1),'BaseValue',mu_factor*Enz_BaseValue,'FaceColor',Enz_color_left1,'EdgeColor','black');
    set(bar1(2),'BaseValue',mu_factor*Enz_BaseValue,'FaceColor',Enz_color_left1_,'EdgeColor','black');
    set(gca,'YScale','log','YColor',Enz_color_left1,'FontSize',20,'FontName','Arial');
    set(gcf,'PaperUnits','centimeters');
    set(gcf,'PaperSize',[20 12]);
    set(gcf,'PaperPosition',[0 0 20 12]);
    
    yyaxis right
    bar2=bar([Enz_right1,Enz_right1_]*mu_factor);
    ylim([mu_factor*Enz_BaseValue,1e-2]);
    ylabel('$E$ ($\mu$mol/m$^2$)','FontSize',20,'FontName','Arial','Interpreter','latex');
    set(bar2(1),'BaseValue',mu_factor*Enz_BaseValue,'FaceColor',Enz_color_right1,'EdgeColor','black');
    set(bar2(2),'BaseValue',mu_factor*Enz_BaseValue,'FaceColor',Enz_color_right1_,'EdgeColor','black');
    set(gca,'YScale','log','YColor',Enz_color_right1,'FontSize',20,'FontName','Arial');
end
if is_mc_feasible==1
    hh9=figure;
    mu_factor=1e0;
    yyaxis left
    bar1=bar([v,v_]*mu_factor);
    ylim([mu_factor*v_BaseValue,1e1]);
    set(gca,'XTick',1:n_rxn,'XTickLabel',rxn_id);
    xtickangle(gca,90);
    ylabel('$v$ (M/s)','FontSize',20,'FontName','Arial','Interpreter','latex');
    set(bar1(1),'BaseValue',mu_factor*v_BaseValue,'FaceColor',v_color_left,'EdgeColor','black');
    set(bar1(2),'BaseValue',mu_factor*v_BaseValue,'FaceColor',v_color_left_,'EdgeColor','black');
    set(gca,'YScale','log','YColor',v_color_left,'FontSize',20,'FontName','Arial');
    set(gcf,'PaperUnits','centimeters');
    set(gcf,'PaperSize',[20 12]);
    set(gcf,'PaperPosition',[0 0 20 12]);
    
    yyaxis right
    bar2=bar([v_right,v_right_]*mu_factor);
    ylim([mu_factor*v_BaseValue,1e-8]);
    ylabel('$v$ (mol/m$^2/s$)','FontSize',20,'FontName','Arial','Interpreter','latex');
    set(bar2(1),'BaseValue',mu_factor*Enz_BaseValue,'FaceColor',v_color_right,'EdgeColor','black');
    set(bar2(2),'BaseValue',mu_factor*Enz_BaseValue,'FaceColor',v_color_right_,'EdgeColor','black');
    set(gca,'YScale','log','YColor',v_color_right,'FontSize',20,'FontName','Arial');
end
if is_sc_feasible==1    
    hh10=figure;
    mu_factor=1e0;
    yyaxis left
    bar1=bar([v1,v1_]*mu_factor);
    ylim([mu_factor*v_BaseValue,1e1]);
    set(gca,'XTick',1:n_rxn,'XTickLabel',rxn_id);
    xtickangle(gca,90);
    ylabel('$v$ (M/s)','FontSize',20,'FontName','Arial','Interpreter','latex');
    set(bar1(1),'BaseValue',mu_factor*v_BaseValue,'FaceColor',v_color_left1,'EdgeColor','black');
    set(bar1(2),'BaseValue',mu_factor*v_BaseValue,'FaceColor',v_color_left1_,'EdgeColor','black');
    set(gca,'YScale','log','YColor',v_color_left1,'FontSize',20,'FontName','Arial');
    set(gcf,'PaperUnits','centimeters');
    set(gcf,'PaperSize',[20 12]);
    set(gcf,'PaperPosition',[0 0 20 12]);
    
    yyaxis right
    bar2=bar([v_right1,v_right1_]*mu_factor);
    ylim([mu_factor*v_BaseValue,1e-7]);
    ylabel('$v$ (mol/m$^2/s$)','FontSize',20,'FontName','Arial','Interpreter','latex');
    set(bar2(1),'BaseValue',mu_factor*Enz_BaseValue,'FaceColor',v_color_right1,'EdgeColor','black');
    set(bar2(2),'BaseValue',mu_factor*Enz_BaseValue,'FaceColor',v_color_right1_,'EdgeColor','black');
    set(gca,'YScale','log','YColor',v_color_right1,'FontSize',20,'FontName','Arial');
end

hh11=figure;
bar1=bar(assim_rate,'stacked');
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
if is_mc_feasible==1 && is_sc_feasible==1
    print(hh1,'-dpdf',strPlotFilePDF1);
    print(hh2,'-dpdf',strPlotFilePDF2);
    print(hh3,'-dpdf',strPlotFilePDF3);
end
if is_mc_km_feasible==1 && is_sc_km_feasible==1
    print(hh4,'-dpdf',strPlotFilePDF4);
    print(hh5,'-dpdf',strPlotFilePDF5);
    print(hh6,'-dpdf',strPlotFilePDF6);
end
if is_mc_feasible==1
    print(hh7,'-dpdf',strPlotFilePDF7);
    print(hh9,'-dpdf',strPlotFilePDF9);
end
if is_sc_feasible==1
    print(hh8,'-dpdf',strPlotFilePDF8);
    print(hh10,'-dpdf',strPlotFilePDF10);
end
print(hh11,'-dpdf',strPlotFilePDF11);
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

function kappa_atps=Fun_ATPS_Sat(x,x_hc,x_hp,atps_par,i_atp,i_adp,i_pi)
kon_pi=atps_par(1);
kon_adp=atps_par(2);
kon_atp=atps_par(3);
kon_hp=atps_par(4);
kon_hc=atps_par(5);
k_dt=atps_par(6);
Kd_pi=atps_par(7);
Kd_adp=atps_par(8);
Kd_atp=atps_par(9);
Kd_hp1=atps_par(10);
Kd_hc1=atps_par(11);
Keq_dt=atps_par(12);
Kd_hp2=atps_par(13);
Kd_hc2=atps_par(14);
Kd_hp3=atps_par(15);
Kd_hc3=atps_par(16);
n_step=atps_par(17);

alpha_p=zeros(n_step,1);
alpha_n=zeros(n_step,1);
alpha_r=zeros(n_step,1);

x_atp=x(i_atp);
x_adp=x(i_adp);
x_pi=x(i_pi);

alpha_p(1)=kon_adp*x_adp;
alpha_p(2)=kon_pi*x_pi;
alpha_p(3)=kon_hp*x_hp;
alpha_p(4)=kon_hc*Kd_hc1;
alpha_p(5)=kon_hp*x_hp;
alpha_p(6)=kon_hc*Kd_hc2;
alpha_p(7)=kon_hp*x_hp;
alpha_p(8)=kon_hc*Kd_hc3;
alpha_p(9)=k_dt;
alpha_p(10)=kon_atp*Kd_atp;

alpha_n(1)=kon_adp*Kd_adp;
alpha_n(2)=kon_pi*Kd_pi;
alpha_n(3)=kon_hp*Kd_hp1;
alpha_n(4)=kon_hc*x_hc;
alpha_n(5)=kon_hp*Kd_hp2;
alpha_n(6)=kon_hc*x_hc;
alpha_n(7)=kon_hp*Kd_hp3;
alpha_n(8)=kon_hc*x_hc;
alpha_n(9)=k_dt/Keq_dt;
alpha_n(10)=kon_atp*x_atp;

alpha_r(1)=alpha_n(n_step)/alpha_p(n_step);
for i=2:n_step
    alpha_r(i)=alpha_n(i-1)/alpha_p(i-1);
end

% Saturation efficiency
sum1=0;
for i=1:n_step
    sum2=0;
    for j=0:n_step-2
        prod=1;
        for k=0:j
            m=i-j+k;
            if m<1
                m=m+n_step;
            end
            prod=prod*alpha_r(m);
        end
        sum2=sum2+prod;
    end
    sum1=sum1+(1+sum2)/alpha_p(i);
end
kappa_atps=1/sum1;
return

function [f,df]=Fun_Obj(ln_x,S,Keq,kapp,varth,rxn,j_crowd,n_crowd,n_met)
Xi=ones(n_crowd,1);
dGamma=ones(n_crowd,n_met);
ddf=zeros(n_crowd,n_met);

x=exp(ln_x);
f=0;
for j=1:n_crowd
    for i=1:n_met
        if S(i,j_crowd(j))~=0
            Xi(j)=Xi(j)*x(i)^S(i,j_crowd(j));
        end
    end
    Xi(j)=Xi(j)/Keq(j_crowd(j));

    f=f+rxn{j_crowd(j)}.MW*varth(j_crowd(j))/kapp(j_crowd(j))/(1-Xi(j));
end

for j=1:n_crowd
    for r=1:n_met
        for i=1:n_met
            if i~=r && S(i,j_crowd(j))~=0
                dGamma(j,r)=dGamma(j,r)*x(i)^S(i,j_crowd(j));
            end
        end
        dGamma(j,r)=dGamma(j,r)*S(r,j_crowd(j))*x(r)^(S(r,j_crowd(j))-1);
    end
end

for j=1:n_crowd
    for r=1:n_met
        ddf(j,r)=rxn{j_crowd(j)}.MW*varth(j_crowd(j))*(dGamma(j,r)/Keq(j_crowd(j)))/(1-Xi(j_crowd(j)))^2/kapp(j_crowd(j));
    end
end
df=sum(ddf,1)';
return

function [f,df]=Fun_Obj_km(ln_x,S,Km,Keq,kapp,varth,rxn,j_crowd,n_crowd,n_met)
Xi=ones(n_crowd,1);
zeta=zeros(n_crowd,1);
prod_s=ones(n_crowd,1);
prod_p=ones(n_crowd,1);
dprod_s=zeros(n_crowd,n_met);
dprod_p=zeros(n_crowd,n_met);
dGamma=ones(n_crowd,n_met);
dzeta=zeros(n_crowd,n_met);
ddf=zeros(n_crowd,n_met);

x=exp(ln_x);
f=0;
for j=1:n_crowd
    for i=1:n_met
        if S(i,j_crowd(j))~=0
            Xi(j)=Xi(j)*x(i)^S(i,j_crowd(j));
        end
    end
    Xi(j)=Xi(j)/Keq(j_crowd(j));

    n_s=0;
    n_p=0;
    for i=1:n_met
        if S(i,j_crowd(j))<0 && Km(i,j_crowd(j))>0
            prod_s(j)=prod_s(j)*(x(i)/Km(i,j_crowd(j)))^(-S(i,j_crowd(j)));
            n_s=n_s+1;
        elseif S(i,j_crowd(j))>0 && Km(i,j_crowd(j))>0
            prod_p(j)=prod_p(j)*(x(i)/Km(i,j_crowd(j)))^S(i,j_crowd(j));
            n_p=n_p+1;
        else
            continue;
        end
    end
    if n_s==0
       prod_s(j)=0;
    end
    if n_p==0
       prod_p(j)=0;
    end
    zeta(j)=prod_s(j)/(1+prod_s(j)+prod_p(j));

    f=f+rxn{j_crowd(j)}.MW*varth(j_crowd(j))/kapp(j_crowd(j))/zeta(j)/(1-Xi(j));
end

for j=1:n_crowd
    for r=1:n_met
        for i=1:n_met
            if i~=r && S(i,j_crowd(j))~=0
                dGamma(j,r)=dGamma(j,r)*x(i)^S(i,j_crowd(j));
            end
        end
        dGamma(j,r)=dGamma(j,r)*S(r,j_crowd(j))*x(r)^(S(r,j_crowd(j))-1);
    end
end

for j=1:n_crowd
    for r=1:n_met
        if S(r,j_crowd(j))<0 && Km(r,j_crowd(j))>0
            dprod_s(j,r)=1;
            for i=1:n_met
                if i~=r && S(i,j_crowd(j))<0 && Km(i,j_crowd(j))>0
                    dprod_s(j,r)=dprod_s(j,r)*(x(i)/Km(i,j_crowd(j)))^(-S(i,j_crowd(j)));
                end
            end
            dprod_s(j,r)=-dprod_s(j,r)*(S(r,j_crowd(j))/Km(r,j_crowd(j)))*(x(r)/Km(r,j_crowd(j)))^(-S(r,j_crowd(j))-1);
        elseif S(r,j_crowd(j))>0 && Km(r,j_crowd(j))>0
            dprod_p(j,r)=1;
            for i=1:n_met
                if i~=r && S(i,j_crowd(j))>0 && Km(i,j_crowd(j))>0
                    dprod_p(j,r)=dprod_p(j,r)*(x(i)/Km(i,j_crowd(j)))^S(i,j_crowd(j));
                end
            end
            dprod_p(j,r)=dprod_p(j,r)*(S(r,j_crowd(j))/Km(r,j_crowd(j)))*(x(r)/Km(r,j_crowd(j)))^(S(r,j_crowd(j))-1);
        else
            continue;
        end
    end
end

for j=1:n_crowd
    for r=1:n_met
        if S(r,j_crowd(j))<0
            dzeta(j,r)=dprod_s(j,r)*(1+prod_p(j))/(1+prod_s(j)+prod_p(j))^2;
        elseif S(r,j_crowd(j))>0
            dzeta(j,r)=-dprod_p(j,r)*prod_s(j)/(1+prod_s(j)+prod_p(j))^2;
        else
            continue;
        end
    end
end

for j=1:n_crowd
    for r=1:n_met
        ddf(j,r)=-rxn{j_crowd(j)}.MW*varth(j_crowd(j))*(dzeta(j,r)*(1-Xi(j_crowd(j)))-zeta(j)*dGamma(j,r)/Keq(j_crowd(j)))/zeta(j)^2/(1-Xi(j_crowd(j)))^2/kapp(j_crowd(j));
    end
end
df=sum(ddf,1)';
return






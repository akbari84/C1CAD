function [is_mc_feasible,uc,xc]=Fun_FeasabilityMC_rGly(pHc,beta,x_c_H2,x_ub,index)
clc;

%% Path strings
strDataPath=Fun_DataPath();
strReadFileMAT1=strcat(strDataPath,'\MAT\parameters-',num2str(index),'.mat');

%% Load data
load(strReadFileMAT1,'met','rxn','n_met','n_rxn','S','Km');

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

kapp=zeros(n_rxn,1);
w=zeros(n_rxn_c,1);
Xi=zeros(n_rxn,1);

%% Parameters
T=25+273.15;                                        % K
pH=[7;pHc;8.0;7.29;7];                              % -             1->e, 2->c, 3->f, 4->single compartment,  5->reference for flux dir correction        
pH_p=2;                                             % 
Is=[0.1;0.1;0.1;0.1;0.1];                           % M             1->e, 2->c, 3->f, 4->single compartment,  5->reference for flux dir correction 
Is_p=0.1;                                           % M 

DPsi=-0.1;                                          % V 
DE_FD=-0.399;                                       % V             DE_FD:=E_FDrd-E_FDox    
DG00w=-238.7;                                       % kJ/mol 
Enz_tot_v=340;                                      % g/L           Based on E. coli cytoplasmatic density 

x_CO2_Q(i_co2)=0.01;                                % M           The concentration of species in the CO2 feed stream (i.e. Q); solubality at 25C = 0.003  
x_H2_f(i_h2)=0.01;                                  % M           The concentration of species in the H2 feed stream (i.e. H2 to f-reactor); solubality at 25C = 0.001 
beta_f=beta/2;                                      % -           beta_f:=FF/c  
beta_e=beta/2;                                      % -           beta_e:=EE/c
eta=0.99;                                           % -           Degree of separation for unit operation at outlet node: fl('CC')*x(pyr,'CC')=eta*fl('c')*x(pyr,'c')  

%% Dependent parameters
DG00_FD=-F*DE_FD;                                   % kJ/mol 
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

%% Log Keq's
ln_Keq=log(Keq(jc));

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
sc=Sc*varth(1:n_rxn_c);
sf=S(:,jf);
se_v=sum(S(:,je_v),2);
se_a=S(:,je_a);

%% Multi-compartment characteristic concentrations
x_mc_scale=(1-beta_f-beta_e)/(3/x_CO2_Q(i_co2)+2*rxn{j_atps}.nH_transport/x_H_p+5/x_H2_f(i_h2));         % Concentration scale of the multi-compartment reactor system
H2_mc=beta_f*x_c_H2/(beta_f+5*x_mc_scale/x_H2_f(i_h2))/(beta_f+beta_e);                                  % H2 concentration in the f stream of the multi-compartment reactor system

%% Reactor optimization subject to crowding constraint without Km (Multi-compartment)
NN=3;
np_rlx=200;                                             % Number of points used to approximat ...<ln(th0-th1/x) by several ...<a_i*ln(x)+b_i constraints linear in ln(x)
d_rlx=5e-8;                                             % Margin from lower/upper bounds of concentration of oxidized cofactors used to generate points for constraint relaxation
eps_rlx=1e-3;                                           % Margin from the relaxed constraint within which solutions are not accepted. This is to avoid solutions where the hydrogenase ractions of the redox compartment are exactly at equilibrium
n_rxn_cf=n_rxn_c+n_rxn_f*np_rlx;
x_lb=1e-7;
x_lb_atp=1e-3;
x_lb_adp=1e-4;
x_lb_pi=1e-3;
x_lb_ppi=1e-3;
%x_ub=0.5;

Aeq=zeros(n_rxn_cf+8,n_met);
beq=zeros(n_rxn_cf+8,1);
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
Aeq(1:n_rxn_c,1:n_met)=Sc';
ii=1;
for i=1:n_rxn_f
    for j=1:np_rlx
        Aeq(n_rxn_c+ii,i_rd(i))=1;
        Aeq(n_rxn_c+ii,i_ox(i))=-1-a_rlx(i,j);
        beq(n_rxn_c+ii)=b_rlx(i,j)-eps_rlx;
        ii=ii+1;
    end
end
Aeq(n_rxn_cf+1,i_nadh)=-1;
Aeq(n_rxn_cf+1,i_nad)=1;
Aeq(n_rxn_cf+2,i_nadph)=-1;
Aeq(n_rxn_cf+2,i_nadp)=1;
Aeq(n_rxn_cf+3,i_Fdrd)=-1;
Aeq(n_rxn_cf+3,i_Fdox)=1;
Aeq(n_rxn_cf+4,i_TRXrd)=-1;
Aeq(n_rxn_cf+4,i_TRXox)=1;
Aeq(n_rxn_cf+5,i_atp)=-1;
Aeq(n_rxn_cf+5,i_adp)=1;
Aeq(n_rxn_cf+5,i_pi)=1;
Aeq(n_rxn_cf+6,i_atp)=-1;
Aeq(n_rxn_cf+6,i_amp)=1;
Aeq(n_rxn_cf+6,i_ppi)=1;
Aeq(n_rxn_cf+7,i_pyr)=1;
Aeq(n_rxn_cf+8,i_h2)=1;
beq(1:n_rxn_c)=ln_Keq;
beq(n_rxn_cf+7)=log(sc(i_pyr)*x_mc_scale/eta);
beq(n_rxn_cf+8)=log(x_c_H2);
Q=a2*Sc*W_*Sc';
c=a1*w_'*Sc'-2*a2*ln_Keq'*W_*Sc';

model.Q=sparse(Q);
model.obj=c;
model.A=sparse(Aeq);
model.rhs=beq;
model.lb=lb1;
model.ub=ub1;
model.modelsense='min';
model.sense=strcat(repmat('<',1,n_rxn_cf),'<<<<<<==');
params.NonConvex=2;
sol=gurobi(model,params);
if ~strcmp(sol.status,'OPTIMAL')
    is_mc_feasible=0;
    uc=NaN;
else 
    is_mc_feasible=1;  
    xc=exp(sol.x);
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
end
return
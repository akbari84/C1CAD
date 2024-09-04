function [is_sc_feasible,uc1,xc1]=Fun_FeasabilitySC_rGly(pHc,beta,x_c_H2,x_ub,index)
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

comp_j=[2;2;2;2;2;2;2;2;2;2;3;3;3;3;1;1;1];
jc=[1;2;3;4;5;6;7;8;9;10];
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

kapp1=zeros(n_rxn,1);
w1=zeros(n_rxn1,1);
Xi1=zeros(n_rxn,1);

%% Parameters
T=25+273.15;                                        % K
pH=[7;6;8.0;pHc;7];                                 % -             1->e, 2->c, 3->f, 4->single compartment,  5->reference for flux dir correction        
pH_p=2;                                             % -
Is=[0.1;0.1;0.1;0.1;0.1];                           % M             1->e, 2->c, 3->f, 4->single compartment,  5->reference for flux dir correction 
Is_p=0.1;                                           % M 

DPsi=-0.1;                                          % V 
DE_FD=-0.399;                                       % V             DE_FD:=E_FDrd-E_FDox    
DG00w=-238.7;                                       % kJ/mol 
Enz_tot_v=340;                                      % g/L           Based on E. coli cytoplasmatic density 

x_CO2_Q(i_co2)=0.01;                                % M             The concentration of species in the CO2 feed stream (i.e. Q); solubality at 25C = 0.003  
x_H2_f(i_h2)=0.01;                                  % M             The concentration of species in the H2 feed stream (i.e. H2 to f-reactor); solubality at 25C = 0.001 
beta_r=beta;                                        % -             beta_r:=RR/c
eta=0.99;                                           % -             Degree of separation for unit operation at outlet node: fl('CC')*x(pyr,'CC')=eta*fl('c')*x(pyr,'c')  

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
for j=1:n_rxn
    kapp1(j)=rxn{j}.kcat;
end

%% Log Keq's
eps=1e-7;
ln_Keq1=[log(Keq1(jc));log(Keq1(jf));log(Keq1(je_v))];
ln_Keq1(j_fdrd)=log(Keq1(j_fdrd)-eps);

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
x_sc_scale=(1-beta_r)/(3/x_CO2_Q(i_co2)+2*rxn{j_atps}.nH_transport/x_H_p+5/x_H2_f(i_h2));                % Concentration scale of the single-compartment reactor system

%% Reactor optimization subject to crowding constraint without Km (Single-compartment)
NN=3;
x_lb=1e-7;
x_lb_atp=1e-3;
x_lb_adp=1e-4;
x_lb_pi=1e-3;
x_lb_ppi=1e-3;

Aeq=zeros(n_rxn1+8,n_met);
beq=zeros(n_rxn1+8,1);
lb1=zeros(n_met,1);
ub1=zeros(n_met,1);

a1=NN*(NN+1)/2;
a2=NN*(2*NN^2+3*NN+1)/12;

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
Q=a2*S1*W1_*S1';
c=a1*w1_'*S1'-2*a2*ln_Keq1'*W1_*S1';
Aeq(1:n_rxn1,:)=S1';
Aeq(n_rxn1+1,i_nadh)=-1;
Aeq(n_rxn1+1,i_nad)=1;
Aeq(n_rxn1+2,i_nadph)=-1;
Aeq(n_rxn1+2,i_nadp)=1;
Aeq(n_rxn1+3,i_Fdrd)=-1;
Aeq(n_rxn1+3,i_Fdox)=1;
Aeq(n_rxn1+4,i_TRXrd)=-1;
Aeq(n_rxn1+4,i_TRXox)=1;
Aeq(n_rxn1+5,i_atp)=-1;
Aeq(n_rxn1+5,i_adp)=1;
Aeq(n_rxn1+5,i_pi)=1;
Aeq(n_rxn1+6,i_atp)=-1;
Aeq(n_rxn1+6,i_amp)=1;
Aeq(n_rxn1+6,i_ppi)=1;
Aeq(n_rxn1+7,i_pyr)=1;
Aeq(n_rxn1+8,i_h2)=1;
beq(1:n_rxn1)=ln_Keq1;
beq(n_rxn1+7)=log(ss1(i_pyr)*x_sc_scale/eta);
beq(n_rxn1+8)=log(x_c_H2);

model.Q=sparse(Q);
model.obj=c;
model.A=sparse(Aeq);
model.modelsense='min';
model.sense=strcat(repmat('<',1,n_rxn1),'<<<<<<==');
model.rhs=beq;
model.lb=lb1;
model.ub=ub1;
params.NonConvex=2;
params.BarHomogeneous=1;
sol=gurobi(model,params);
if strcmp(sol.status,'OPTIMAL')
    is_sc_feasible=1;
    xc1=exp(sol.x);
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
else 
    is_sc_feasible=0;  
    uc1=NaN;
    xc1=zeros(n_met,1);
end
return
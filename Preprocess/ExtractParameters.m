function ExtractParameters()
clc;

%% Parameters to set up path strings 
index=1;

%% Path strings
strDataPath=Fun_DataPath();
strReadFileEXL=strcat(strDataPath,'\EXL\Network-',num2str(index),'.xlsx');
strSaveFileMAT=strcat(strDataPath,'\MAT\parameters-',num2str(index),'.mat');

%% Read Excel files
met_kegg=readcell(strReadFileEXL,'Sheet','met','Range','B:B');
met_id=readcell(strReadFileEXL,'Sheet','met','Range','C:C');
met_name=readcell(strReadFileEXL,'Sheet','met','Range','D:D');
pKa_str=readcell(strReadFileEXL,'Sheet','met','Range','F:F');
z_A=readmatrix(strReadFileEXL,'Sheet','met','Range','G:G');
NH_A=readmatrix(strReadFileEXL,'Sheet','met','Range','H:H');
DfG00_A=readmatrix(strReadFileEXL,'Sheet','met','Range','I:I');
is_aq=readmatrix(strReadFileEXL,'Sheet','met','Range','J:J');

rxn_kegg=readcell(strReadFileEXL,'Sheet','rxn','Range','B:B');
rxn_id=readcell(strReadFileEXL,'Sheet','rxn','Range','C:C');
rxn_gene=readcell(strReadFileEXL,'Sheet','rxn','Range','D:D');
rxn_ec=readcell(strReadFileEXL,'Sheet','rxn','Range','E:E');
rxn_name=readcell(strReadFileEXL,'Sheet','rxn','Range','F:F');
MW=readmatrix(strReadFileEXL,'Sheet','rxn','Range','G:G');
S_w=readmatrix(strReadFileEXL,'Sheet','rxn','Range','H:H');
nH_transport=readmatrix(strReadFileEXL,'Sheet','rxn','Range','I:I');
kcat=readmatrix(strReadFileEXL,'Sheet','rxn','Range','J:J');
DrG0tr_ex=readmatrix(strReadFileEXL,'Sheet','rxn','Range','K:K');
dir=readmatrix(strReadFileEXL,'Sheet','rxn','Range','L:L');
dir_kin=readmatrix(strReadFileEXL,'Sheet','rxn','Range','M:M');
is_mem=readmatrix(strReadFileEXL,'Sheet','rxn','Range','N:N');

S=readmatrix(strReadFileEXL,'Sheet','S');
Km=readmatrix(strReadFileEXL,'Sheet','Km');

%% Adjusting the dimension of S and Km
[n_met,n_rxn]=size(S);
S=S(2:n_met,2:n_rxn);
Km=Km(2:n_met,2:n_rxn);
met_id=met_id(2:n_met);
rxn_id=rxn_id(2:n_rxn);
n_met=n_met-1;
n_rxn=n_rxn-1;

%% Create met object
met=cell(n_met,1);
for i=1:n_met
    met{i}.kegg=met_kegg{i+1};
    met{i}.id=met_id{i};                % first row already removed 
    met{i}.name=met_name{i+1};
    met{i}.pKa=str2num(pKa_str{i+1});
    met{i}.z_A=z_A(i+1);
    met{i}.NH_A=NH_A(i+1);
    met{i}.DfG00_A=DfG00_A(i+1);
    met{i}.N=size(met{i}.pKa,2);
    met{i}.is_aq=is_aq(i+1);            % is in aqueous phase?
end

%% Create rxn object
rxn=cell(n_rxn,1);
for i=1:n_rxn
    rxn{i}.kegg=rxn_kegg{i+1};
    rxn{i}.id=rxn_id{i};                % first row already removed 
    rxn{i}.gene=rxn_gene{i+1};
    rxn{i}.ec=rxn_ec{i+1};
    rxn{i}.name=rxn_name{i+1};
    rxn{i}.MW=MW(i+1);
    rxn{i}.S_w=S_w(i+1);
    rxn{i}.nH_transport=nH_transport(i+1);
    rxn{i}.kcat=kcat(i+1);
    rxn{i}.DrG0tr_ex=DrG0tr_ex(i+1);
    rxn{i}.dir=dir(i+1);
    rxn{i}.dir_kin=dir_kin(i+1);
    rxn{i}.is_mem=is_mem(i+1);          % is membrane protein?
end

%% Save MAT
save(strSaveFileMAT,'met','rxn','n_met','n_rxn','S','Km','met_id','rxn_id');
return
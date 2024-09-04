function CalcFormationEnergy()
clc;

%% Parameters to set up path strings 
index=101;

%% Path strings
strDataPath=Fun_DataPath();
strResultsPathEXL=Fun_ExcelPath();
strReadFileMAT1=strcat(strDataPath,'\MAT\parameters-',num2str(index),'.mat');
strSaveFileMAT=strcat(strDataPath,'\MAT\formation_energy-',num2str(index),'.mat');
strSaveFileEXL=strcat(strResultsPathEXL,'\formation_energy-',num2str(index),'.xlsx');

%% Load MAT files
load(strReadFileMAT1,'met','n_met','met_id');

%% Constants
Rgas=8.31446261815324;                              % J/K/mol 

%% Parameters  
T=25+273.15;                                        % K

%% Memory allocation
n00=zeros(n_met,1);
n00_max=0;
for i=1:n_met
    n00(i)=size(met{i}.pKa,2)+1; %number of charge states (H bound, Mg bound, and the respective DG00)
    if n00(i)>n00_max
        n00_max=n00(i);
    end
end
DG00a=zeros(n_met,n00_max);   %H bound
Ka=zeros(n_met,n00_max);  
KKa=zeros(n_met,n00_max);   

%% Calculate formation energies of charged states
Ka(:,1)=1;
KKa(:,1)=1;
for i=1:n_met
    DG00a(i,1)=met{i}.DfG00_A;
    for k=2:n00(i)
        DG00a(i,k)=DG00a(i,k-1)-0.001*Rgas*T*log(10)*met{i}.pKa(k-1);
    end
    for k=2:n00(i)
        Ka(i,k)=10^(-met{i}.pKa(k-1));
    end
    for k=2:n00(i)
        KKa(i,k)=KKa(i,k-1)*Ka(i,k);
    end
end

%% Write to Excel
%sheet1
tempcell=convertCharsToStrings('met_id');
xlRange=strcat(xlscol(1),'1');
writematrix(tempcell,strSaveFileEXL,'sheet','H','range',xlRange);
xlRange=strcat('A2:A',num2str(n_met+1));
writecell(met_id,strSaveFileEXL,'sheet','H','range',xlRange);

for k=1:n00_max
    tempcell=convertCharsToStrings(strcat('DG00_',num2str(k-1)));
    xlRange=strcat(xlscol(k+1),'1');
    writematrix(tempcell,strSaveFileEXL,'sheet','H','range',xlRange);
    xlRange=strcat(xlscol(k+1),'2:',xlscol(k+1),num2str(n_met+1));
    writematrix(DG00a(:,k),strSaveFileEXL,'sheet','H','range',xlRange);
end

%sheet2
tempcell=convertCharsToStrings('met_id');
xlRange=strcat(xlscol(1),'1');
writematrix(tempcell,strSaveFileEXL,'sheet','Ka','range',xlRange);
xlRange=strcat('A2:A',num2str(n_met+1));
writecell(met_id,strSaveFileEXL,'sheet','Ka','range',xlRange);

for k=1:n00_max
    tempcell=convertCharsToStrings(strcat('Ka_',num2str(k-1)));
    xlRange=strcat(xlscol(k+1),'1');
    writematrix(tempcell,strSaveFileEXL,'sheet','Ka','range',xlRange);
    xlRange=strcat(xlscol(k+1),'2:',xlscol(k+1),num2str(n_met+1));
    writematrix(Ka(:,k),strSaveFileEXL,'sheet','Ka','range',xlRange);
end

%sheet3
tempcell=convertCharsToStrings('met_id');
xlRange=strcat(xlscol(1),'1');
writematrix(tempcell,strSaveFileEXL,'sheet','KKa','range',xlRange);
xlRange=strcat('A2:A',num2str(n_met+1));
writecell(met_id,strSaveFileEXL,'sheet','KKa','range',xlRange);

for k=1:n00_max
    tempcell=convertCharsToStrings(strcat('KKa_',num2str(k-1)));
    xlRange=strcat(xlscol(k+1),'1');
    writematrix(tempcell,strSaveFileEXL,'sheet','KKa','range',xlRange);
    xlRange=strcat(xlscol(k+1),'2:',xlscol(k+1),num2str(n_met+1));
    writematrix(KKa(:,k),strSaveFileEXL,'sheet','KKa','range',xlRange);
end

%sheet8
tempcell=convertCharsToStrings('met_id');
xlRange=strcat(xlscol(1),'1');
writematrix(tempcell,strSaveFileEXL,'sheet','n00','range',xlRange);
xlRange=strcat('A2:A',num2str(n_met+1));
writecell(met_id,strSaveFileEXL,'sheet','n00','range',xlRange);

tempcell=convertCharsToStrings('n_states');
xlRange=strcat(xlscol(2),'1');
writematrix(tempcell,strSaveFileEXL,'sheet','n00','range',xlRange);
xlRange=strcat(xlscol(2),'2:',xlscol(2),num2str(n_met+1));
writematrix(n00(:),strSaveFileEXL,'sheet','n00','range',xlRange);

%% Save
save(strSaveFileMAT,'DG00a','n00','n00_max');

return
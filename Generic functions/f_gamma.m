function f=f_gamma(Ic_s,A_PHI_PZ_s,B_PHI_PZ_s)
    %f=-A_PHI_PZ_s*(sqrt(Ic_s)/(1+B_PHI_PZ_s*sqrt(Ic_s))+2*log(1+B_PHI_PZ_s*sqrt(Ic_s))/B_PHI_PZ_s); %Pitzar
    f=-3*A_PHI_PZ_s*sqrt(Ic_s)/(1+(4/3)*B_PHI_PZ_s*sqrt(Ic_s)); %Extended Debye-Huckel
return
function f=f_gamma_EDH(Ic,A_PHI_PZ,B_PHI_PZ)
    f=3*A_PHI_PZ*sqrt(Ic)/(1+(4/3)*B_PHI_PZ*sqrt(Ic)); %Extended Debye-Huckel
return
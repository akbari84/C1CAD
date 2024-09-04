function f=f_phi(Ic_s,A_PHI_PZ_s,B_PHI_PZ_s)
    f=-A_PHI_PZ_s*sqrt(Ic_s)/(1+B_PHI_PZ_s*sqrt(Ic_s));
return
function [dtruss_dA1, dtruss_dA2, dtruss_dA3, dtruss_dA4, dtruss_dA5] = truss_compliance_sensitivity_by_area_by_set(Kee_a1,Kee_a2,Kee_a3,Kee_a4,Kee_a5, Fe, B_a1,B_a2,B_a3,B_a4,B_a5, u_c)
dtruss_dA1 = (Fe'*((Kee_a1\Fe) - B_a1*u_c))/10e11;
dtruss_dA2 = (Fe'*((Kee_a2\Fe) - B_a2*u_c))/10e11;
dtruss_dA3 = (Fe'*((Kee_a3\Fe) - B_a3*u_c))/10e11;
dtruss_dA4 = (Fe'*((Kee_a4\Fe) - B_a4*u_c))/10e11;
dtruss_dA5 = (Fe'*((Kee_a5\Fe) - B_a5*u_c))/10e11;
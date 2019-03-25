function [dtruss_compliance_sensitivity_dA] = truss_compliance_sensitivity(Kee_a, Fe, B_a, u_c)
dtruss_compliance_sensitivity_dA = (Fe'*((Kee_a\Fe) - B_a*u_c))/10e15;
a=5;


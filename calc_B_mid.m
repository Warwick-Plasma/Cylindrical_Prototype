% calc_B_mid.m
% Average the values of the new magnetic fields, and the one from the
% previous step. The fields calculated here will be evaluated at the same
% time as the eletric field

Bxm_mid = 0.5 * (Bxm + Bxm_old);
Brm_mid = 0.5 * (Brm + Brm_old);
Btm_mid = 0.5 * (Btm + Btm_old);
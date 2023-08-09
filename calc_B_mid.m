% finalise_fields.m
% Average the values of the new magnetic fields, and the one from the
% previous step. The fields calculated here will be evaluated at the same
% time as the eletric field. Also clear the current densities ahead of the
% next particle push

% Get fields at the electric field evaluation time
Bxm_mid = 0.5 * (Bxm + Bxm_old);
Brm_mid = 0.5 * (Brm + Brm_old);
Btm_mid = 0.5 * (Btm + Btm_old);

% Reset current densities
Jxm = zeros(nx+6, nr+5, n_mode);
Jrm = zeros(nx+5, nr+6, n_mode);
Jtm = zeros(nx+5, nr+5, n_mode);
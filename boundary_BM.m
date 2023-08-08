% Field stagger in the top-r ghost-cell
%
%   x x x x x x x x x x x x x x x x x x x x x x x x x                             
%  Bxm(i,nr+6,:)                                    |
%  Jrm(i,nr+6,:)       Btm(i+1,nr+6,:)              |
%  Erm(i,nr+6,:)                                    |
%   |                                               |
%   |                  Brm(i+1,nr+5,:)              |
%  Jtm(i,nr+5,:)-------Jxm(i+1,nr+5,:)--------------|----  r = (nr+2)*dr
%  Etm(i,nr+5,:)       Exm(i+1,nr+5,:)              |
%   |                                               |
%  Bxm(i,nr+5,:)                                    |
%  Jrm(i,nr+5,:)       Btm(i+1,nr+5,:)              |
%  Erm(i,nr+5,:)                                    |
%   |                                               |
%
% Equations on this boundary come from the Faraday-Lenz equation,
% evaluating gradients at the position and time of Etm(i,nr+5) (for Bxm)
% and Exm(i,nr+5) (for Btm). The conditions Ex=-Bt and Et=Bx have also been
% used. Here, (Kx,Kr)=(0,1) and (1,0) correspond to two different
% approaches to solving the same equations, which the code may interpolate
% between. The Smilie default (0,1) introduces terms from the
% Ampere-Maxwell equation

% Vary Kx and Kr to get different methods to solve the same equation
Kx = 0;
Kr = 1;
cosphi = Kr / sqrt(Kx^2 + Kr^2);
CB_BM  = cosphi/(1.0 + cosphi); % Theta is always taken equal to zero.
CE_BM  = 1 - CB_BM;

% Boundary is on the highest primal r point
one_ov_rlocal = 1/(dr*(nr+2));

% Coeffs for Bxm
factor= 1 / (1 + dt/dr);
Alpha_Bl_Rmax = (1 - dt/dr) * factor;
Beta_Bl_Rmax = CE_BM * dt * one_ov_rlocal * factor;
Gamma_Bl_Rmax = CB_BM * dt / dx *factor;

% Coeffs for Btm
factor = 1 / (1 + dt/dr + 0.5 * CB_BM * dt * one_ov_rlocal);
Alpha_Bt_Rmax = (-1 + dt/dr - 0.5 * CB_BM * dt * one_ov_rlocal) * factor;
Beta_Bt_Rmax = (1 - dt/dr - 0.5 * CB_BM * dt * one_ov_rlocal) * factor;
Gamma_Bt_Rmax = (1 + dt/dr - 0.5 * CB_BM * dt * one_ov_rlocal) * factor;
Epsilon_Bt_Rmax = dt * one_ov_rlocal * factor;
Delta_Bt_Rmax = dt/dx * factor;

for imode=0:Nmode
              
    % For Bxm^(p,d)
    % SMILEI ERROR: MISSING +2*CB_BM*DT*FACTOR_BL*JXM(I,NR+5,IMODE+1)
    for i=1:nx+5
        Bxm(i,nr+6,imode+1) = Bxm_old(i,nr+5,imode+1) ...
            - Alpha_Bl_Rmax * (Bxm(i,nr+5,imode+1) - Bxm_old(i,nr+6,imode+1)) ...
            + Gamma_Bl_Rmax * (Brm(i+1,nr+5,imode+1) + Brm_old(i+1,nr+5,imode+1) - Brm(i,nr+5,imode+1) - Brm_old(i,nr+5,imode+1)) ...
            - Beta_Bl_Rmax * imagi * imode * (Erm(i,nr+6,imode+1) + Erm(i,nr+5,imode+1)) ...
            - 2 * Beta_Bl_Rmax * Etm(i,nr+5,imode+1);
    end
    
    % For Btm^(d,d)
    % SMILEI ERROR: MISSING -2*CB_BM*DT*FACTOR_BT*JTM(I,NR+5,IMODE+1)
    % SMILEI ERROR: SHOULD BE A PLUS IN (Brm(i,nr+5,imode+1) + Brm_old(i,nr+5,imode+1))
    for i=2:nx+6
        Btm(i,nr+6,imode+1) = Alpha_Bt_Rmax * Btm(i,nr+5,imode+1) ...
            + Beta_Bt_Rmax  * Btm_old(i,nr+6,imode+1) ...
            + Gamma_Bt_Rmax * Btm_old(i,nr+5,imode+1) ...
            - imagi * imode * CB_BM * Epsilon_Bt_Rmax  * (Brm(i,nr+5,imode+1) - Brm_old(i,nr+5,imode+1)) ...
            - CE_BM * Delta_Bt_Rmax * (Erm(i,nr+6,imode+1) + Erm(i,nr+5,imode+1) - Erm(i-1,nr+6,imode+1) - Erm(i-1,nr+5,imode+1)) ;
    end
end
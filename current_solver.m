% void ProjectorAM2Order::currents(   ElectroMagnAM *emAM, 
%                                     Particles &particles, 
%                                     unsigned int ipart,
%                                     double invgf, 
%                                     int *iold, 
%                                     double *deltaold, 
%                                     std::complex<double> *array_eitheta_old, 
%                                     bool diag_flag, int ispec)

% Current stagger in the bottom-left cell in the simulation grid
%   |                            |
%   |                            |
%   |----------------------------|----
%   |                            |
%  Jrm(3,4,:)                    |
%   |                            |
%   |                            |
%  Jtm(3,3,:)---Jxm(4,3,:)-------|----

% Macro-particle incides for currents
% Vertical bars:    cell edges
% Points:           Jx or Jr evaluation points
% Cross:            particle position before push
% Horizontal lines: primal cells, numbered
% Roman numerals:   indices of Jl_p, and Jr_p
% Jr          I     II   III    IV    V
% Jx    I     II   III    IV    V
%    |     |     |   x |     |     |     |
%    |  .  |  .  |  .  |  .  |  .  |  .  |  -> x or r
%         ___   ___   ___   ___   ___
%          1     2     3     4     5
% Macro-particle always starts in primal cell 3, macro-particle shape can 
% only pass 4 J evaluation points. Left of primal 1 and right of primal 5 
% are un-reachable in one timestep. This is why Jl_p(1) and Jr_p(5) are
% always 0.

%% Variable declaration & initialization

% Multiply this by 1/r to get charge*weight/cell_volume for a cell centred
% at r
QWr_over_V = charge(ipart) * weight(ipart) / (2*pi*dx*dr);

% Collect constant terms for calculating current densities
crl_p = QWr_over_V * dx / dt;
crr_p = QWe_over_V / dt;

% Arrays used for the Esirkepov projection method
Sl0 = zeros(1,5);
Sl1 = zeros(1,5);
Sr0 = zeros(1,5);
Sr1 = zeros(1,5);
DSl = zeros(1,5);
DSr = zeros(1,5);
Jl_p = zeros(1,5);
Jr_p = zeros(1,5);
C_m = 1;

%% Locate particles & calculate Esirkepov coef. S, DS and W

% Locate the particle on the primal grid at former time-step & calculate coeff. S0
% The terms deltax and deltar were previously calculated in
% fields_to_particles.m, called before particle_push.m
delta = deltax;
delta2 = delta^2;
Sl0(2) = 0.5 * (delta2 - delta + 0.25);
Sl0(3) = 0.75 - delta2;
Sl0(4) = 0.5 * (delta2 + delta + 0.25);

delta = deltar;
delta2 = delta^2;
Sr0(2) = 0.5 * (delta2 - delta + 0.25);
Sr0(3) = 0.75 - delta2;
Sr0(4) = 0.5 * (delta2 + delta + 0.25);

% Calculate exponential coefficients

% exp(i*theta_old), where exp_m_theta = exp(-i*theta_old), from
% fields_to_particles
theta_old = conj(exp_m_theta); 

% exp(i*theta_new)
r_new = sqrt(pos_y(ipart)^2 + pos_z(ipart)^2);
eitheta = (pos_y(ipart) + Icpx * pos_z(ipart)) / r_new;

% exp(i * (theta_new - theta_old)/2)
% std::sqrt and MATLAB sqrt keeps the root with positive real part which is what we need here.
e_delta_m1 = sqrt(eitheta * conj(theta_old)); 

% exp(i * (theta_new + theta_old)/2)
e_bar_m1 = theta_old * e_delta_m1;

% Locate the particle on the primal grid at current time-step & calculate coeff. S1
xpn = pos_x(ipart) / dx;
ipo = ip;
ip = round(xpn);
ip_m_ipo = ip - ipo;
delta = xpn - ip;
delta2 = delta^2;
Sl1(ip_m_ipo+2) = 0.5 * (delta2 - delta + 0.25);
Sl1(ip_m_ipo+3) = 0.75 - delta2;
Sl1(ip_m_ipo+4) = 0.5 * (delta2 + delta + 0.25);

ypn = r_new / dr;
jpo = jp;
jp = round(ypn);
jp_m_jpo = jp - jpo;
delta  = ypn - jp;
delta2 = delta^2;
Sr1(jp_m_jpo+2) = 0.5 * (delta2 - delta + 0.25);
Sr1(jp_m_jpo+3) = 0.75 - delta2;
Sr1(jp_m_jpo+4) = 0.5 * (delta2 + delta + 0.25);

% Calculate differences in cell-weights over the time-step. Here, Sl0(2)
% and Sr0(2) refer to the primal cell (centred on the bottom-left corner of
% the normal grid), which contains the macro-particle at the old position
for i = 1:5
    DSl(i) = Sl1(i) - Sl0(i);
    DSr(i) = Sr1(i) - Sr0(i);
end

% r at t = t0 - dt/2
% I THINK SMILEI HAS THIS LINE WRONG, AND I THINK IT SHOULD ACTUALLY BE
% r_bar = ((jpo + deltar)*dr + r_new) * 0.5;
% Where (jpo + deltar)*dr is equivalent to r_old
r_bar = (jpo*dr + deltar + r_new) * 0.5; 

% 1 / r corresponding to the radius of the centres of the primal cells
% (normal cell edges). invR_local(3) is 1/r for the primal-centre closest
% to the old r-position
cell_r = abs((jpo-2)*dr:dr:(jpo+2)*dr);
for ir = 1:5
    if cell_r(ir) == 0
        % Smilei gives comment: "No Verboncoeur correction"
        invR_local(ir) = 8.0 / dr;
    else
        invR_local(ir) = 1.0 / cell_r(ir);
    end
end

% 1 / r corresponding to the radius of primal cell-edges, shifted +0.5*dr
% from invR_local
invRd = 1.0 / (cell_r + 0.5*dr);

%% Pre-calculate terms for current calculation

% Initial value of crt_p for imode = 0.
% Q * W * r * p_theta / gamma / V
crt_p = QWr_over_V * (pz(ipart)*real(e_bar_m1) - py(ipart)*imag(e_bar_m1)) * inv_gf(ipart);

% Compute everything which has no theta dependence
tmpJl = zeros(1,5);
for j = 1:5
    % Q * W * [mean weight-fraction in r] / (2 * pi * r * dr * dt)
    % Approximately: Q * W * [mean-weight-fraction-in-r] / A / dt
    tmpJl(j) = crl_p * (Sr0(j) + 0.5*DSr(j)) * invR_local(j);
end

% Weight change past Jx evaluation points. See header diagram for
% explanation why Jl_p(1) is always 0
Jl_p(1)= 0.;
for i = 2:5
    % Jl_p(1) is zero from the geometry of the system
    % Any change in cell 2 must be from current passing Jl_p(2)
    % If no change in 3, an equal current flows through Jl_p(3) into Jl_p(2)
    Jl_p(i)= Jl_p(i-1) - DSl(i-1);
end

% Vd(i): for Jr evaluation point i, gives ratio between radii at Jr(i+1), Jr(i)
% tmpJr(i): Q * W / (2*pi*r*dx) / dt * [r-weight-diff-between Jr(i) and Jr(i+1)]
% tmpJr(5) = 0 the same reason Jl_p(1) = 0, no weight can pass that point
Vd = zeros(1,5);
tmpJr = zeros(1,5);
for j=-4:-1
    Vd(-j) = (cell_r(-j) + dr) * invRd(-j);
    tmpJr(-j) = crr_p * DSr(-j+1) * invRd(-j) * dr;
end

% Current past each r evaluation point, neglecting x-weight-fraction
Jr_p(5) = 0.;
for j=-4:-1
    % - Jr_p(5) is zero from the geometry of the system
    % - Any change in cell 5 is from current passing Jr_p(4)
    % - If no weight change in cell 4, then equal current flows between
    %   Jr_p(3) and Jr_p(4) BUT area at III is smaller than IV by factor 
    %   r3/r4, so current density greater by r4/r3 (this is Vd). tmpJr 
    %   gives current due to excess weight-change.
    Jr_p(-j) =  Jr_p(-j+1) * Vd(-j) + tmpJr(-j);
end

% Set initial e_delta and e_delta_inv values for m=0 special case Jtheta 
% calculation. These will be overriden with correct values after m=0
% I THINK SMILEI HAS THIS LINE WRONG, AND I THINK IT SHOULD ACTUALLY BE
% e_delta = 0.5;
e_delta = 1.5;
e_delta_inv = 0.5;

% exp(i * m * (theta_new + theta_old)/2) for m=0
e_bar = 1.0;

% Compute division by R in advance for Jt and rho evaluation.
for j = 1:5
    Sr0(j) = Sr0(j) * invR_local(j);
    Sr1(j) = Sr1(j) * invR_local(j);
end

%% Perform current calculation

% The current modes calculated in this section are missing a 1/(2*pi) 
% prefactor 

for imode = 0:Nmode-1
    % Method used for m=0 Jt is different to m>0 Jt, here we set values for
    % m>0
    if imode > 0
        e_delta = e_delta * e_delta_m1;
        e_bar = e_bar * e_bar_m1;
        e_delta_inv = 1.0 / e_delta - 1.0;
        
        % New Jt pre-factor - calculation for m>0 assumes small dtheta
        crt_p = charge_weight * Icpx * e_bar / (dt * imode) * 2.0 * r_bar;

        % Multiply modes > 0 by 2 and C_m = 1 otherwise.
        C_m = 2.0 * e_bar; 
    end

    % Add contribution J_p to global array
    % Jl^(d,p)
    for i = 2:5
        for j=1:5
            % C_m (Lifschitz) * [x-weight-change] * [perp. current density]
            Jxm(ipo+i,jpo+j,imode+1) = Jxm(ipo+i,jpo+j,imode+1) + C_m * Jl_p(i) * tmpJl(j);
        end
    end

    % Jr^(p,d)
    for i = 1:5
        for j = 2:5
            % C_m (Lifschitz) * [mean-x-weight] * [r current density]
            Jrm(ipo+i,jpo+j,imode+1) = Jrm(ipo+i,jpo+j,imode+1) + C_m * (Sl0(i) + 0.5*DSl(i)) * Jr_p(j) ;
        end 
    end

    % Jt^(p,p)
    for i = 1:5
        for j = 1:5
            % When m=0, [crt_p: Q * W * r * p_theta / gamma / V] * [0.5 * (cell_weight_new + cell_weight_old) / r] => only if typo is present
            % When m>0, small angle approx reduces expression to: Cm * [QW/V] * [v_theta: mean_r * dtheta / dt] * [0.5 * (cell_weight_new + cell_weight_old)]
            Jtm(ipo+i,jpo+j,imode+1) = Jtm(ipo+i,jpo+j,imode+1) + crt_p * (Sr1(j)*Sl1(j)*e_delta_inv - Sr0(j)*Sl0(i)*(e_delta-1.0));
        end
    end

    % Restore e_delta correct initial value, after m=0 special case
    if imode == 0 
        e_delta = 1.0;
    end
end
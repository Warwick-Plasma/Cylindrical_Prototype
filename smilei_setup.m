% smilei_setup.m
% Create field arrays, load particles

imagi = sqrt(-1);

% Initialise field and current mode arrays
Jxm = zeros(nx+6, nr+5, n_mode);
Jrm = zeros(nx+5, nr+6, n_mode);
Jtm = zeros(nx+5, nr+5, n_mode);

Exm = zeros(nx+6, nr+5, n_mode);
Erm = zeros(nx+5, nr+6, n_mode);
Etm = zeros(nx+5, nr+5, n_mode);

Bxm = zeros(nx+5, nr+6, n_mode);
Brm = zeros(nx+6, nr+5, n_mode);
Btm = zeros(nx+6, nr+6, n_mode);

Bxm_old = zeros(nx+5, nr+6, n_mode);
Brm_old = zeros(nx+6, nr+5, n_mode);
Btm_old = zeros(nx+6, nr+6, n_mode);

Bxm_mid = zeros(nx+5, nr+6, n_mode);
Brm_mid = zeros(nx+6, nr+5, n_mode);
Btm_mid = zeros(nx+6, nr+6, n_mode);

% Geometric grid parameters
cell_r_min = 0:dr:(nr-1)*dr;
cell_r_max = dr:dr:nr*dr;
cell_vol = pi*dx*(cell_r_max.^2 - cell_r_min.^2);

% Laser parameters
Amp_y_xmin = cos(pol);
Amp_z_xmin = sin(pol);
Amp_y_xmax = cos(pol);
Amp_z_xmax = sin(pol);
t_centre = 0.5 * t_fwhm; 

% Calculate particle count
ix_min = max(floor(el_x_min / dx) + 1, 1);
ix_max = min(floor(el_x_max / dx) + 1, nx);
ir_min = max(floor(el_r_min / dr) + 1, 1);
ir_max = min(floor(el_r_max / dr) + 1, nr);
npart = (ix_max - ix_min + 1) * (ir_max - ir_min + 1) * el_ppc;

% Initialise particle arrays
px = zeros(1,npart);
py = zeros(1,npart);
pz = zeros(1,npart);
pos_x = zeros(1,npart);
pos_y = zeros(1,npart);
pos_z = zeros(1,npart);
weight = zeros(1,npart);
invgf = zeros(1,npart);   % 1/gamma (Lorentz)
in_sim = true(1,npart);

% Fields on particles
Ex = zeros(1,npart);
Ey = zeros(1,npart);
Ez = zeros(1,npart);
Bx = zeros(1,npart);
By = zeros(1,npart);
Bz = zeros(1,npart);

% All particles are electrons, charge -1, mass 1 in SMILEI units
charge = -1;
mass = 1;

% Set positions of macro-particles
ipart = 1;
for ix = ix_min:ix_max
    for ir = ir_min:ir_max
        for i = 1:el_ppc
            pos_x(ipart) = (ix-rand())*dx;

            pos_r = (ir-rand())*dr;
            pos_th = 2*pi*rand();
            pos_y(ipart) = pos_r * cos(pos_th);
            pos_z(ipart) = pos_r * sin(pos_th);

            % Let particle weight scale linearly with r
            weight(ipart) = pos_r;
            ipart = ipart+1;
        end

        % Normalise weights in current cell to sum to expected real
        % particle count
        real_part_no = el_ne * cell_vol(ir);
        norm_factor = real_part_no / sum(weight(ipart-el_ppc:ipart-1));
        weight(ipart-el_ppc:ipart-1) = ...
            weight(ipart-el_ppc:ipart-1) * norm_factor;
    end
end

% Set momenta of macro-particles
fluct_p_mag = el_temp * rand(1,npart);
fluct_theta = acos(1 - 2*rand(1,npart));
fluct_phi = 2*pi*rand(1,npart);
gamma_drift = sqrt(1/(1-el_drift_vx^2));
px_drift = gamma_drift * el_drift_vx;
px = px_drift + fluct_p_mag.*sin(fluct_theta).*cos(fluct_phi);
py = fluct_p_mag.*sin(fluct_theta).*sin(fluct_phi);
pz = fluct_p_mag.*cos(fluct_theta);

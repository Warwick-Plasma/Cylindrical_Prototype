% run_prototype.m
%
% The main script which takes the user specification, and runs 
% the cylindrical PIC prototype. This code is designed to model an electron
% bunch passing through a Gaussian laser beam
%
% This code uses the normalised units of SMILEI

%% User input

% Control block
dx = 0.05;          % Cell size in x
dr = 0.05;          % Cell size in r
nx = 200;           % Number of cells in the x direction
nr = 100;           % Number of cells in the r direction
dt = 0.015;         % Simulation time-step
t_end = 15;         % Run-time of simulation
n_mode = 2;         % Number of axial modes to model

% Laser block (all lasers share same parameters)
las_x_min = true;            % There is a Gaussian beam laser on x_min
las_x_max = false;           % There is a Gaussian beam laser on x_max
omega = 6.28;                % Angular frequency
waist = 5;                   % r-coord where field drops to 1/e at focus
foc_x = 5;                   % Distance from boundary to focal spot in x
foc_r = 0;                   % r-position of focal spot
t_fwhm = 12;                 % FWHM of temporal Gaussian profile   
a0 = 0.00854;                % Laser pump strength
pol = 0;                     % Polarisation: 0 for Ey, pi/2 for Ez

% Electron species block
el_ne = 0.000353;           % Electron number density (constant)
el_x_min = 9.9;             % Initial low-x boundary of e- bunch
el_x_max = 10.0;            % Initial high-x boundary of e- bunch
el_r_min = 0.0;             % Initial low-r boundary of e- bunch
el_r_max = 0.1;             % Initial high-r boundary of e- bunch
el_ppc = 5;                 % Number of electrons per cell
el_drift_vx = -0.999975;    % Drift of electrons [units of c]
el_temp = 2.0e-3;           % Random p-drift between 0 and this [unit me*c]

%% Set-up

smilei_setup;

%% PIC loop

t = 0;
while t < t_end

    t = t + dt

    % Update particle positions and currents
    for ipart = 1:npart
        if (~in_sim(ipart))
            % Ignore particles outside the simulation window
            continue
        end
        fields_to_particles;
        particle_push;
        current_solver;
        current_BCs;
    end

    % Update fields
    save_B_old;
    E_update;
    B_update;
    boundary_SM;
    boundary_BM;
    calc_B_mid;

end
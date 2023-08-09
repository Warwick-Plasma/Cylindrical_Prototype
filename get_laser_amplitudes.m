% get_laser_amplitudes.m
% Calculates the By and Bz values of the incident laser pulse, whether on
% the x_min or x_max boundary. We assume all lasers are Gaussian beams
% which share the same distance to focus

% Derived parameters (amplitudes)
if x_min_boundary
    if las_x_min
        amplitudeY = Amp_y_xmin * a0 * omega;
        amplitudeZ = Amp_z_xmax * a0 * omega;
    else
        % No laser on this boundary
        las_By = 0;
        las_Bz = 0;
        return
    end
else
    if las_x_max
        amplitudeY = Amp_y_xmin * a0 * omega; 
        amplitudeZ = Amp_y_xmax * a0 * omega;
    else
        % No laser on this boundary
        las_By = 0;
        las_Bz = 0;
        return
    end
end

% Derived parameters (Gaussian beam)
Zr = omega * waist^2 * 0.5;
w  = sqrt(1 / (1 + (foc_x/Zr)^2));
invWaist2 = (w/waist)^2;
coeff = -omega * foc_x * w^2 / (2 * Zr^2);
t_sigma = (0.5 * t_fwhm)^2 / log(2);

% Calculate profiles
r_env = w * exp(-invWaist2 * (rj - foc_r)^2);
t_env = exp(-(t - t_centre)^2 / t_sigma);
phase = coeff * (rj - foc_r)^2;

las_By = amplitudeZ * r_env * t_env * sin(omega*t - phase);
las_Bz = amplitudeY * r_env * t_env * sin(omega*t - phase);
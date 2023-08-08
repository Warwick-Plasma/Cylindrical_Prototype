% Derived parameters (amplitudes)
if las_x_min
    amplitudeY = Amp_y_xmin * a0 * omega;
    amplitudeZ = Amp_z_xmax * a0 * omega;
else
    amplitudeY = Amp_y_xmin * a0 * omega; 
    amplitudeZ = Amp_y_xmax * a0 * omega;
end

% Derived parameters (Gaussian beam)
Zr = omega * waist^2 * 0.5;
w  = sqrt(1 / (1 + (foc_x/Zr)^2));
invWaist2 = (w/waist)^2;
coeff = -omega * foc_x * w^2 / (2 * Zr^2);
t_sigma = (0.5 * t_fwhm)^2 / log(2);

% Calculate profiles
r_env = w * exp(-invWaist2 * (r - foc_r)^2);
t_env = exp(-(t - t_centre)^2 / t_sigma);
phase = coeff * (r - foc_r)^2;

las_By = amplitudeZ * r_env * t_env * sin(omega*t - phase);
las_Bz = amplitudeY * r_env * t_env * sin(omega*t - phase);
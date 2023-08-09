% el_x_px_phase_space.m
% Plot a heatmap of the electron phase space in x and px.
% This script assumes you have already run a cylindrical simulation

px_bin_count = 50;
x_bin_count = nx;

% Set phase-space limits
x_edges = linspace(0,nx*dx,x_bin_count+1);
px_edges = sign(px_drift) * linspace(0.5*abs(px_drift), 1.5*abs(px_drift), px_bin_count + 1);
if sign(px_drift) == -1
    px_edges = flip(px_edges);
end

% Find phase space bin indices of each particle
npart_plot = length(in_sim(in_sim==1));
[~,~,~,bin_x,bin_px] = histcounts2(pos_x(in_sim),px(in_sim),x_edges,px_edges);

% Count macro-particle weight in each bin
phase_space = zeros(x_bin_count, px_bin_count);
for i = 1:npart_plot 
    phase_space(bin_x(i),bin_px(i)) = phase_space(bin_x(i),bin_px(i)) + weight(i);
end

% Plot heatmap
x_centres = 0.5*(x_edges(2:end) + x_edges(1:end-1));
px_centres = 0.5*(px_edges(2:end) + px_edges(1:end-1));
[x_plot, px_plot] = meshgrid(x_centres, px_centres);
surf(x_plot, px_plot,phase_space','EdgeColor','none');
cbar = colorbar;
view(2);
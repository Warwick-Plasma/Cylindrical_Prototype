% plot_Er.m
% Calculates Er along a specified theta from the available field modes, and
% creates a heatmap plot

theta = 0;

Er = zeros(nx,nr);
for im = 0:n_mode-1
    Er = Er + real(Erm(4:nx+3,4:nr+3,im+1))*cos(im*theta) + imag(Erm(4:nx+3,4:nr+3,im+1))*sin(im*theta);
end

x_edges = linspace(0,nx*dx,nx+1);
r_edges = linspace(0,nr*dr,nr+1);
x_centres = 0.5*(x_edges(2:end) + x_edges(1:end-1));
r_centres = 0.5*(r_edges(2:end) + r_edges(1:end-1));
[x_plot, r_plot] = meshgrid(x_centres, r_centres);
surf(x_plot, r_plot, Er', 'EdgeColor', 'none');
cbar = colorbar;
view(2);
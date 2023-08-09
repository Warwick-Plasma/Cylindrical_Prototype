% particle_push.m
% This script is called by the PIC loop, along with a particle index ipart.
% Here, we use the interpolated fields from fields_to_particles.m to 
% perform a Boris particle-push, updating both position and momentum for 
% particle ipart. We also tag particles leaving the simulation window for
% later removal

charge_over_mass_dts2 = charge * 0.5 * dt / mass;

% Initiate half-acceleration in the electric field
pxsm = charge_over_mass_dts2 * Ex(ipart);
pysm = charge_over_mass_dts2 * Ey(ipart);
pzsm = charge_over_mass_dts2 * Ez(ipart);

umx = px(ipart) + pxsm;
umy = py(ipart) + pysm;
umz = pz(ipart) + pzsm;

% Rotation in the magnetic field
local_invgf = charge_over_mass_dts2 / sqrt( 1.0 + umx*umx + umy*umy + umz*umz );
Tx = local_invgf * Bx(ipart);
Ty = local_invgf * By(ipart);
Tz = local_invgf * Bz(ipart);
inv_det_T = 1.0/(1.0 + Tx*Tx + Ty*Ty + Tz*Tz);

pxsm = pxsm + ((1.0+Tx*Tx-Ty*Ty-Tz*Tz)*umx + 2.0*(Tx*Ty+Tz)*umy + 2.0*(Tz*Tx-Ty)*umz) * inv_det_T;
pysm = pysm + (2.0*(Tx*Ty-Tz)*umx + (1.0-Tx*Tx+Ty*Ty-Tz*Tz)* umy  + 2.0*(Ty*Tz+Tx)*umz) * inv_det_T;
pzsm = pzsm + (2.0*(Tz*Tx+Ty)*umx + 2.0*(Ty*Tz-Tx)* umy  + (1.0-Tx*Tx-Ty*Ty+Tz*Tz)*umz) * inv_det_T;

% Finalise half-acceleration in the electric field
local_invgf = 1.0/sqrt(1.0 + pxsm*pxsm + pysm*pysm + pzsm*pzsm);
invgf(ipart) = local_invgf;

px(ipart) = pxsm;
py(ipart) = pysm;
pz(ipart) = pzsm;

% Move the particle
local_invgf = local_invgf * dt;
pos_x(ipart) = pos_x(ipart) + pxsm*local_invgf;
pos_y(ipart) = pos_y(ipart) + pysm*local_invgf;
pos_z(ipart) = pos_z(ipart) + pzsm*local_invgf;

% Mark escaped particles
if pos_x(ipart) > nx*dx || pos_x(ipart) < 0
    in_sim(ipart) = false;
elseif pos_y(ipart)^2 + pos_z(ipart)^2 > (nr*dr)^2 
    in_sim(ipart) = false;
end
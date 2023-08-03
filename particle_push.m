ipart2 = ipart - ipart_buffer_offset;

charge_over_mass_dts2 = charge(ipart) * one_over_mass * dts2;

% Initiate half-acceleration in the electric field
pxsm = charge_over_mass_dts2 * Ex(ipart2);
pysm = charge_over_mass_dts2 * Ey(ipart2);
pzsm = charge_over_mass_dts2 * Ez(ipart2);

umx = momentum_x(ipart) + pxsm;
umy = momentum_y(ipart) + pysm;
umz = momentum_z(ipart) + pzsm;

% Rotation in the magnetic field
local_invgf = charge_over_mass_dts2 / sqrt( 1.0 + umx*umx + umy*umy + umz*umz );
Tx = local_invgf * Bx(ipart2);
Ty = local_invgf * By(ipart2);
Tz = local_invgf * Bz(ipart2);
inv_det_T = 1.0/(1.0 + Tx*Tx + Ty*Ty + Tz*Tz);

pxsm = pxsm + ((1.0+Tx*Tx-Ty*Ty-Tz*Tz)*umx + 2.0*(Tx*Ty+Tz)*umy + 2.0*(Tz*Tx-Ty)*umz) * inv_det_T;
pysm = pysm + (2.0*(Tx*Ty-Tz)*umx + (1.0-Tx*Tx+Ty*Ty-Tz*Tz)* umy  + 2.0*(Ty*Tz+Tx)*umz) * inv_det_T;
pzsm = pzsm + (2.0*(Tz*Tx+Ty)*umx + 2.0*(Ty*Tz-Tx)* umy  + (1.0-Tx*Tx-Ty*Ty+Tz*Tz)*umz) * inv_det_T;

% Finalise half-acceleration in the electric field
local_invgf = 1.0/sqrt(1.0 + pxsm*pxsm + pysm*pysm + pzsm*pzsm);
invgf(ipart2) = local_invgf;

momentum_x(ipart) = pxsm;
momentum_y(ipart) = pysm;
momentum_z(ipart) = pzsm;

% Move the particle
local_invgf = local_invgf * dt;
position_x(ipart) = position_x(ipart) + pxsm*local_invgf;
position_y(ipart) = position_y(ipart) + pysm*local_invgf;
position_z(ipart) = position_z(ipart) + pzsm*local_invgf;

% Mark escaped particles
if position_x(ipart) > x_max || position_x(ipart) < x_min
    in_sim(ipart) = false;
elseif position_y(ipart)^2 + position_z(ipart)^2 > r_max^2 
    in_sim(ipart) = false;
end
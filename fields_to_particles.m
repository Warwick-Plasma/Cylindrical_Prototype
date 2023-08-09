% fields_to_particles.m
% This script is called by the PIC loop, along with a particle index ipart.
% Here, we interpolate the electric and magnetic fields to the particle
% position, to get Ex,Er and Et (labelled Ex, Ey, Ez). Then we perform a
% rotation to obtain Ex, Ey and Ez.

% Treat mode 0 first

% Normalized particle position
xpn = pos_x(ipart) / dx;
r = sqrt(pos_y(ipart)^2 + pos_z(ipart)^2);
rpn = r / dr;
exp_m_theta = (pos_y(ipart) - imagi*pos_z(ipart)) / r;  % exp(-i*theta)
exp_mm_theta = 1;                                       % exp(-i*m*theta)

% Indices of the central cells on primal and staggered grids
% Primal grid: cell between -0.5*dx to 0.5*dx has index 3
% Dual grid: cell between 0 and dx has index 4
ip = round(xpn) + 3;
id = round(xpn+0.5) + 3;
jp = round(rpn) + 3;
jd = round(rpn+0.5) + 3;

% Calculation of the coefficient for interpolation
deltax = xpn - (id-3) + 0.5;
delta2 = deltax^2;
coeff_xd(1) = 0.5 * (delta2 - deltax + 0.25);
coeff_xd(2) = 0.75 - delta2;
coeff_xd(3) = 0.5 * (delta2 + deltax + 0.25);

deltax = xpn - (ip-3);
delta2 = deltax^2;
coeff_xp(1) = 0.5 * (delta2 - deltax + 0.25);
coeff_xp(2) = 0.75 - delta2;
coeff_xp(3) = 0.5 * (delta2 + deltax + 0.25);

deltar = rpn - (jd-3) + 0.5;
delta2  = deltar^2;
coeff_yd(1) = 0.5 * (delta2 - deltar + 0.25);
coeff_yd(2) = 0.75 - delta2;
coeff_yd(3) = 0.5 * (delta2 + deltar + 0.25);

deltar = rpn - (jp-3);
delta2 = deltar^2;
coeff_yp(1) = 0.5 * (delta2 - deltar + 0.25);
coeff_yp(2) = 0.75 - delta2;
coeff_yp(3) = 0.5 * (delta2 + deltar + 0.25);
   
% Here we assume that mode 0 is real
Ex(ipart) = 0.0;
Ey(ipart) = 0.0;
Ez(ipart) = 0.0;
Bx(ipart) = 0.0;
By(ipart) = 0.0;
Bz(ipart) = 0.0;
for idx = -1:1
    for idy = -1:1
        % Interpolation of Exm(d,p), Erm(p,d), Etm(p,p)
        Ex(ipart) = Ex(ipart) + coeff_xd(2+idx) * coeff_yp(2+idy) * Exm(id+idx,jp+idy,1);
        Ey(ipart) = Ey(ipart) + coeff_xp(2+idx) * coeff_yd(2+idy) * Erm(ip+idx,jd+idy,1);
        Ez(ipart) = Ez(ipart) + coeff_xp(2+idx) * coeff_yp(2+idy) * Etm(ip+idx,jp+idy,1);
        % Interpolation of Bxm(p,d), Brm(d,p), Btm(d,d)
        Bx(ipart) = Bx(ipart) + coeff_xp(2+idx) * coeff_yd(2+idy) * Bxm_mid(ip+idx,jd+idy,1);
        By(ipart) = By(ipart) + coeff_xd(2+idx) * coeff_yp(2+idy) * Brm_mid(id+idx,jp+idy,1);
        Bz(ipart) = Bz(ipart) + coeff_xd(2+idx) * coeff_yd(2+idy) * Btm_mid(id+idx,jd+idy,1);
    end
end

% Sum real parts of higher field modes
if n_mode >= 2
for im = 2:n_mode      
    % Generate successive exp(i*m*theta) by multiplying by exp(i*theta)
    exp_mm_theta = exp_mm_theta * exp_m_theta;
    
    Ex_sum = 0;
    Ey_sum = 0;
    Ez_sum = 0;
    Bx_sum = 0;
    By_sum = 0;
    Bz_sum = 0;
    for idx = -1:1
        for idy = -1:1
            % Interpolation of Exm(d,p), Erm(p,d), Etm(p,p)
            Ex_sum = Ex_sum + coeff_xd(2+idx) * coeff_yp(2+idy) * Exm(id+idx,jp+idy,im);
            Ey_sum = Ey_sum + coeff_xp(2+idx) * coeff_yd(2+idy) * Erm(ip+idx,jd+idy,im);
            Ez_sum = Ez_sum + coeff_xp(2+idx) * coeff_yp(2+idy) * Etm(ip+idx,jp+idy,im);
            % Interpolation of Bxm(p,d), Brm(d,p), Btm(d,d)
            Bx_sum = Bx_sum + coeff_xp(2+idx) * coeff_yd(2+idy) * Bxm_mid(ip+idx,jd+idy,im);
            By_sum = By_sum + coeff_xd(2+idx) * coeff_yp(2+idy) * Brm_mid(id+idx,jp+idy,im);
            Bz_sum = Bz_sum + coeff_xd(2+idx) * coeff_yd(2+idy) * Btm_mid(id+idx,jd+idy,im);
        end
    end

    Ex(ipart) = Ex(ipart) + real(Ex_sum * exp_mm_theta);
    Ey(ipart) = Ey(ipart) + real(Ey_sum * exp_mm_theta);
    Ez(ipart) = Ez(ipart) + real(Ez_sum * exp_mm_theta);
    Bx(ipart) = Bx(ipart) + real(Bx_sum * exp_mm_theta);
    By(ipart) = By(ipart) + real(By_sum * exp_mm_theta);
    Bz(ipart) = Bz(ipart) + real(Bz_sum * exp_mm_theta);
end
end
 
% Rotate field into the cartesian y,z coordinates
delta2 = real(exp_m_theta) * Ey(ipart) + imag(exp_m_theta) * Ez(ipart);
Ez(ipart) = -imag(exp_m_theta)*Ey(ipart) + real(exp_m_theta)*Ez(ipart);
Ey(ipart) = delta2;
delta2 = real(exp_m_theta) * By(ipart) + imag(exp_m_theta) * Bz(ipart);
Bz(ipart) = -imag(exp_m_theta)*By(ipart) + real(exp_m_theta)*Bz(ipart);
By(ipart) = delta2;

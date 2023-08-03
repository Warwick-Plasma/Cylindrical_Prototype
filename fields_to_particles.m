Icpx = sqrt(-1);

% Treat mode 0 first

% Normalized particle position
xpn = pos_x(ipart) / dx;
r = sqrt(pos_y(ipart)^2 + pos_z(ipart)^2);
rpn = r / dr;
exp_m_theta = (pos_y(ipart) - Icpx*pos_z(ipart)) / r;  % exp(-i*theta)
exp_mm_theta = 1.;                                     % exp(-i*m*theta)

% Indices of the central cells on primal and staggered grids
ip = round(xpn);
id = round(xpn+0.5);
jp = round(rpn);
jd = round(rpn+0.5);

% Calculation of the coefficient for interpolation
deltax = xpn - id + 0.5;
delta2 = deltax^2;
coeff_xd(1) = 0.5 * (delta2 - deltax + 0.25);
coeff_xd(2) = 0.75 - delta2;
coeff_xd(3) = 0.5 * (delta2 + deltax + 0.25);

deltax = xpn - ip;
delta2 = deltax^2;
coeff_xp(1) = 0.5 * (delta2 - deltax + 0.25);
coeff_xp(2) = 0.75 - delta2;
coeff_xp(3) = 0.5 * (delta2 + deltax + 0.25);

deltar = rpn - jd + 0.5;
delta2  = deltar^2;
coeff_yd(1) = 0.5 * (delta2 - deltar + 0.25);
coeff_yd(2) = 0.75 - delta2;
coeff_yd(3) = 0.5 * (delta2 + deltar + 0.25);

deltar = rpn - jp;
delta2 = deltar^2;
coeff_yp(1) = 0.5 * (delta2 - deltar + 0.25);
coeff_yp(2) = 0.75 - delta2;
coeff_yp(3) = 0.5 * (delta2 + deltar + 0.25);

% First index for summation
ip = ip - i_domain_begin;
id = id - i_domain_begin;
jp = jp - j_domain_begin;
jd = jd - j_domain_begin;
   
% Here we assume that mode 0 is real
Ex(ipart2) = 0.0;
Ey(ipart2) = 0.0;
Ez(ipart2) = 0.0;
Bx(ipart2) = 0.0;
By(ipart2) = 0.0;
Bz(ipart2) = 0.0;
for idx = -1:1
    for idy = -1:1
        % Interpolation of El^(d,p), Er^(p,d), Et^(p,p)
        Ex(ipart2) = Ex(ipart2) + coeff_xd(2+idx) * coeff_yp(2+idy) * El(id+idx,jp+idy,1);
        Ey(ipart2) = Ey(ipart2) + coeff_xp(2+idx) * coeff_yd(2+idy) * Er(id+idx,jp+idy,1);
        Ez(ipart2) = Ez(ipart2) + coeff_xp(2+idx) * coeff_yp(2+idy) * Et(id+idx,jp+idy,1);
        % Interpolation of Bl^(p,d), Br^(d,p), Bt^(d,d)
        Bx(ipart2) = Bx(ipart2) + coeff_xp(2+idx) * coeff_yd(2+idy) * Bl(id+idx,jp+idy,1);
        By(ipart2) = By(ipart2) + coeff_xd(2+idx) * coeff_yp(2+idy) * Br(id+idx,jp+idy,1);
        Bz(ipart2) = Bz(ipart2) + coeff_xd(2+idx) * coeff_yd(2+idy) * Bt(id+idx,jp+idy,1);
    end
end
    
for im = 1:nmode      
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
            Ex_sum = Ex_sum + coeff_xd(2+idx) * coeff_yp(2+idy) * El(id+idx,jp+idy,im+1);
            Ey_sum = Ey_sum + coeff_xp(2+idx) * coeff_yd(2+idy) * Er(id+idx,jp+idy,im+1);
            Ez_sum = Ez_sum + coeff_xp(2+idx) * coeff_yp(2+idy) * Et(id+idx,jp+idy,im+1);
            Bx_sum = Bx_sum + coeff_xp(2+idx) * coeff_yd(2+idy) * Bl(id+idx,jp+idy,im+1);
            By_sum = By_sum + coeff_xd(2+idx) * coeff_yp(2+idy) * Br(id+idx,jp+idy,im+1);
            Bz_sum = Bz_sum + coeff_xd(2+idx) * coeff_yd(2+idy) * Bt(id+idx,jp+idy,im+1);
        end
    end

    Ex(ipart2) = Ex(ipart2) + real(Ex_sum * exp_mm_theta);
    Ey(ipart2) = Ey(ipart2) + real(Ey_sum * exp_mm_theta);
    Ez(ipart2) = Ez(ipart2) + real(Ez_sum * exp_mm_theta);
    Bx(ipart2) = Bx(ipart2) + real(Bx_sum * exp_mm_theta);
    By(ipart2) = By(ipart2) + real(By_sum * exp_mm_theta);
    Bz(ipart2) = Bz(ipart2) + real(Bz_sum * exp_mm_theta);
end
 
% Translate field into the cartesian y,z coordinates
delta2 = real(exp_m_theta) * Ey(ipart2) + imag(exp_m_theta) * Ez(ipart2);
Ez(ipart2) = -imag(exp_m_theta)*Ey(ipart2) + real(exp_m_theta)*Ez(ipart2);
Ey(ipart2) = delta2;
delta2 = real(exp_m_theta) * By(ipart2) + imag(exp_m_theta) * Bz(ipart2);
Bz(ipart2) = -imag(exp_m_theta)*By(ipart2) + real(exp_m_theta)*Bz(ipart2);
By(ipart2) = delta2;

% boundary_SM.m
% Solve boundaries on x_min and x_max using SMILEI's "Silver-Muller"
% method. 
% 
% We calculate the magnetic field using the cylindrical form of the
% Ampere-Maxwell equation, with gradients in time and space centred on
% % Etm(1,j,im+1)    (for Brm(1,j,im+1))
% % Erm(1,j,im+1)    (for Btm(1,j,im+1))
% % Etm(nx+5,j,im+1) (for Brm(nx+6,j,im+1))
% % Erm(nx+5,j,im+1) (for Btm(nx+6,j,im+1))
%
% Field stagger in the bottom-left cell in the simulation grid
%   |                            |
%   |                            |
%   |----------------------------|----
%   |                            |
%   |                            |
%  Bxm(3,4,:)                    |
%  Jrm(3,4,:)   Btm(4,4,:)       |
%  Erm(3,4,:)                    |
%   |                            |
%   |           Brm(4,3,:)       |
%  Jtm(3,3,:)---Jxm(4,3,:)-------|----
%  Etm(3,3,:)   Exm(4,3,:)
%
% For inwards propagating S and P polarised lasers on x_min/max, we have:
% % S = 0.5*(Er + c*Bt)  (x_min)
% % P = 0.5*(Et - c*Br)  (x_min)
% % S = 0.5*(Er - c*Bt)  (x_max)
% % P = 0.5*(Er + c*Br)  (x_max)
%
% As we are dealing with plane waves, these can be rewritten as
% % S = +c*Bt  (x_min)
% % P = -c*Br  (x_min)
% % S = -c*Bt  (x_max)
% % P = +c*Br  (x_max)
%
% From considering Br and Bt expansions, and cylindrical to cartesian
% conversions:
% % Br1 = By + i*Bz
% % Bt1 = Bz - i*By

factor = 1 / (1 + dt/dx);
Alpha_Xmin = 2 * factor;
Beta_Xmin = -(1 - dt/dx) * factor;
Gamma_Xmin = -4 * factor;
Delta_Xmin = -dt/dr * factor;
Epsilon_Xmin = imagi * dt * factor;

% POSSIBLE SMILEI ERROR: I THINK (DELTA_XMAX = +DT/DR * FACTOR) - SEE BELOW
Alpha_Xmax = 2 * factor;
Beta_Xmax = -(1.0 - dt/dx) * factor;
Gamma_Xmax = 4 * factor;
Delta_Xmax = -dt/dr * factor;
Epsilon_Xmax = - imagi * dt * factor;

for im=0:n_mode-1

    %% x_min boundary
    x_min_boundary = true;

    % For Brm^(d,p)
    for j=4:nr+5

        % A single laser  
        laser_P_xmin = 0.0;
        rj = (j-3)*dr;
        if im==1
            get_laser_amplitudes;
            laser_P_xmin = laser_P_xmin - las_By - imagi * las_Bz;
        end
        
        % POSSIBLE SMILEI ERRORS: 
        % 1. THERE IS A SEMI-COLON AFTER BYW, CUTTING IT OFF FROM THE 
        %   DELTA_XMIN LINE. I BELIEVE THIS IS A MISTAKE. I HAVE CORRECTED 
        %   THIS
        % 2. MISSING -DT*JTM(1,J,im+1)/(1+DT/DX)
        % 3. WE DO NOT CONSIDER J=3, BRM ON THE R=0 AXIS
        Brm(1,j,im+1) = Alpha_Xmin * Etm(1,j,im+1) ...
            + Beta_Xmin * Brm(2,j,im+1) ...
            - Gamma_Xmin * laser_P_xmin ...
            + Delta_Xmin * (Bxm(1,j+1,im+1) - Bxm(1,j,im+1));
    end
    
    
    % For Btm^(d,d)
    for j=4:nr+6
             
        % A single laser
        laser_S_xmin = 0;
        rj = (j-3.5)*dr;
        if im==1
            get_laser_amplitudes;
            laser_S_xmin = laser_S_xmin + las_Bz - imagi * las_By;
        end

        % x=Xmin
        Btm(1,j,im+1) = -Alpha_Xmin * Erm(1,j,im+1) ...
            + Beta_Xmin * Btm(2,j,im+1) ...
            + Gamma_Xmin * laser_S_xmin ...
            + Epsilon_Xmin * im / rj * Bxm(1,j,im+1);
    end
    
    % Redo condition on axis for Bt because it was modified
    % Btm(1,3,2): Use d(B_perp)/dtheta = 0 on r=0
    % Btm(1,3,im+1): Use Btm(3)=-2iBrm(3) - Btm(4) where Brm(3)=0
    if im == 1 
        Btm(1,3,2) = -2 * imagi * Brm(1,3,2) - Btm(1,4,2);
    else
        Btm(1,3,im+1) = -Btm(1,4,im+1);
    end

    %% x_max boundary
    x_min_boundary = false;

    % For Brm^(d,p)
    for j=4:nr+5
    
        % A single laser
        laser_P_xmax = 0;
        rj = (j-3)*dr;
        if im==1
            get_laser_amplitudes;
            laser_P_xmax = laser_P_xmax + las_By + imagi * las_Bz;
        end

        % POSSIBLE SMILEI ERRORS: 
        % 1. SHOULDN'T DELTA_XMAX BE +DT/DR INSTEAD OF -DT/DR?
        % 2. MISSING +DT*JTM(1,J,im+1)/(1+DT/DX)
        % 3. WE DO NOT CONSIDER J=3, BRM ON THE R=0 AXIS
        Brm(nx+6,j,im+1) = -Alpha_Xmax * Etm(nx+5,j,im+1) ...
            + Beta_Xmax * Brm(nx+5,j,im+1) ...
            + Gamma_Xmax * laser_P_xmax ...
            + Delta_Xmax * (Bxm(nx+5,j+1,im+1) - Bxm(nx+5,j,im+1));
    end
    
    % For Bt^(d,d)
    for j=4:nr+6
    
        % A single laser
        laser_S_xmax = 0;
        rj = (j-3.5)*dr;
        if im==1
            get_laser_amplitudes;
            laser_S_xmax = laser_S_xmax - las_Bz + imagi * las_By;
        end

        Btm(nx+6,j,im+1) = Alpha_Xmax * Erm(nx+5,j,im+1) ...
            + Beta_Xmax  * Btm(nx+5,j,im+1) ...
            - Gamma_Xmax * laser_S_xmax ...
            + Epsilon_Xmax * im / rj * Bxm(nx+5,j,im+1);
    end

    % Redo condition on axis for Bt because it was modified
    % Btm(nx+6,3,2): Use d(B_perp)/dtheta = 0 on r=0
    % Btm(nx+6,3,im+1): Use Btm(3)=-2iBrm(3) - Btm(4) where Brm(3)=0
    if im == 1
        Btm(nx+6,3,2)= -2 * imagi * Brm(nx+6,3,2) - Btm(nx+6,4,2);
    else
        Btm(nx+6,3,im+1) = -Btm(nx+6,4,im+1);
    end
end
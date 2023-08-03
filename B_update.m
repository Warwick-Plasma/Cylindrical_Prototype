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
%  If we have nx cells, then there are nx dual (cell-centre) points in the
%  grid, and 3 dual points either side in ghost cells (nx+6). There are
%  nx+1 primal points on cell-edges in the window, with 2 more primal
%  points either side (nx+5)

for imode=0:Nmode

    % Magnetic field Bx^(p,d)
    for i=1:nx+5
        for j=4:j<nr+5
            Bxm(i,j,imode+1) = Bxm(i,j,imode+1) - dt/((j-3.5)*dr) ...
                * ((j-3)*Etm(i,j,imode+1) - (j-4)*Etm(i,j-1,imode+1) ...
                + imagi * imode * Erm(i,j,imode+1));
        end
    end
    
    % Magnetic field Br^(d,p)
    for i=2:nx+5
        for j=4:nr+5
            Brm(i,j,imode+1) = Brm(i,j,imode+1) ...
                + dt/dx * (Etm(i,j,imode+1) - Etm(i-1,j,imode+1)) ...
                + imagi * dt * imode / ((j-3)*dr) * Exm(i,j,imode+1);
        end
    end

    % Magnetic field Bt^(d,d)
    for i=2:nx+5
        for j=4:nr+5
            Btm(i,j,imode+1) = Btm(i,j,imode+1) ...
                + dt/dr * (Exm(i,j,imode+1) - Exm(i,j-1,imode+1)) ...
                - dt/dx * (Erm(i,j,imode+1) - Erm(i-1,j,imode+1));
        end
    end
end
    
% Conditions on axis
% m=0
% Brm(:,3,1): Multivalue forbidden
% Brm(:,2,1): Physical field reflection
% Btm(:,3,1): Physical field reflection
% Bxm(:,3,1): Physical field reflection
Brm(:,3,1) = 0;
Brm(:,2,1) = -Brm(:,4,1);
Btm(:,3,1) = -Btm(:,4,1);
Bxm(:,3,1) = Bxm(:,4,1);

% m=1
% Bxm(:,3,2): Physical field reflection
% Brm(2:nx+5,3,2): Solve normally, use one-sided Taylor treatment on Ex/r
% Brm(:,2,2): Physical field reflection
% Btm(:,3,2): Use d(B_perp)/dtheta = 0 on r=0
Bxm(:,3,2) = -Bxm(:,4,2);
Brm(2:nx+5,3,2) = Brm(2:nx+5,3,2) + imagi * dt/dr * Exm(2:nx+5,4,2) ...
    + dt/dx * (Etm(2:nx+5,3,2) - Etm(1:nx+4,3,2));
Brm(:,2,2) = Brm(:,4,2);
Btm(:,3,2) = -2 * imagi * Brm(:,3,2) - Btm(:,4,2);

% m>1
% I BELIEVE THERE IS A MISTAKE IN BRM(2) - WOULD BE BETTER IF PHYSICAL
% REFLECTION
% Bxm(:,3,imode+1): Force-zero at r=0 - multivalue forbidden
% Brm(:,3,imode+1): Multivalue forbidden
% Brm(:,2,imode+1): Physical reflection, correct for even m
% Btm(:,3,imode+1): Use Btm(3)=-2iBrm(3) - Btm(4) where Brm(3)=0
for imode = 2:nmode-1
    Bxm(:,3,imode+1) = -Bxm(:,4,imode+1);
    Brm(:,3,imode+1) = 0;
    Brm(:,2,imode+1) = -Brm(:,4,imode+1);
    Btm(:,3,imode+1) = -Btm(:,4,imode+1);
end
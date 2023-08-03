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
    
    % Electric field Exm(d,p)
    for i=1:nx+6
        for j=4:nr+5
            Exm(i,j,imode+1) = Exm(i,j,imode+1) - dt*Jxm(i,j,imode+1 ) ...
                + dt/((j-3)*dr) * ((j-2.5)*Btm(i,j+1,imode+1) - (j-3.5)*Btm(i,j,imode+1)) ...
                + imagi*dt*imode / ((j-3)*dr) * Brm(i,j,imode+1);
        end
    end

    % Electric field Erm(p,d)
    for i=1:nx+5
        for j=4:nr+6
            Erm(i,j,imode+1) = Erm(i,j,imode+1) - dt*Jrm(i,j,imode+1) ...
                - dt/dx * (Btm(i+1,j,imode+1) - Btm(i,j,imode+1)) ...
                - imagi*dt*imode / ((j-2.5)*dr) * Bxm(i,j,imode+1);                               
        end
    end

    % Electric field Etm(p,p)
    for i=1:nx+5
        for j=4:nr+5
            Etm(i,j,imode+1) = Etm(i,j,imode+1) - dt*Jtm(i,j,imode+1) ...
                + dt/dx * (Brm(i+1,j,imode+1) - Brm(i,j,imode+1)) ...
                - dt/dr * (Bxm(i,j+1,imode+1) - Bxm(i,j,imode+1));
        end
    end
end

% Conditions on axis
% m=0
% Etm(:,3,1): Multivalue forbidden
% Etm(:,2,1): Physical field reflection
% Erm(:,3,1): Physical field reflection
% Exm(:,3,1): Solve normally, 1-sided Taylor treatment on (1/r)*d/dr(r*Bt)
% Exm(:,2,1): Physical field reflection
Etm(:,3,1) = 0;
Etm(:,2,1) = - Etm(:,4,1);
Erm(:,3,1) = - Erm(:,4,1);
Exm(:,3,1) = Exm(:,3,1) + 4*dt/dr*Btm(:,4,1) - dt*Jxm(:,3,1);
Exm(:,2,1) = Exm(:,4,1);

% m=1
% Exm(:,3,2): Multivalue forbidden
% Exm(:,2,2): Physical field reflection
% Etm(:,3,2): From div(E)=0 (rho(m=1)=0), and d(E_perp)/dtheta = 0 for r=0
% Etm(:,2,2): Physical field reflection   
% Erm(:,3,2): Implied from Etm condition
Exm(:,3,2) = 0;
Exm(:,2,2) = - Exm(:,3,2);
Etm(:,3,2) = -imagi/8*(9*Erm(:,4,2) - Erm(:,5,2));
Etm(:,2,2) = Etm(:,4,2);    
Erm(:,3,2)= 2*imagi*Etm(:,3,2) - Erm(:,4,2);

% m > 1
% I BELIEVE SMILEI IS WRONG ON EXM(2),ETM(4). FOR EVEN M, EX(2) = EX(4), 
% AND FOR ODD M, ET(2) = ET(4), AND THE CURRENT SIGNS ARE WRONG
% Exm(:,3,imode+1): Multivalue forbidden
% Exm(:,2,imode+1): Physical field reflection (only valid for odd m)
% Erm(:,4,imode+1): div(E)=0 (rho(m>1)=0), d(E_perp)/dtheta=0, Etm(3)=0
% Erm(:,3,imode+1): Force Er(m>1) = 0, multivalue forbidden
% Etm(:,3,imode+1): Multivalue forbidden
% Etm(:,2,imode+1): Physical field reflection (only valid for even m)
for imode = 2:nmode-1
    Exm(:,3,imode+1) = 0; 
    Exm(:,2,imode+1) = -Exm(:,4,imode+1); 
    Erm(:,4,imode+1) = Erm(:,5,imode+1) / 9.;
    Erm(:,3,imode+1) = -Erm(:,4,imode+1);
    Etm(:,3,imode+1) = 0;
    Etm(:,2,imode+1) = -Etm(:,4,imode+1);
end
% void ProjectorAM2Order::apply_axisBC(std::complex<double> *rhoj,std::complex<double> *Jl, std::complex<double> *Jr, std::complex<double> *Jt, unsigned int imode+1, bool diag_flag )

% Current stagger in the bottom-left cell in the simulation grid
%   |                            |
%   |                            |
%   |----------------------------|----
%   |                            |
%  Jrm(3,4,:)                    |
%   |                            |
%   |                            |
%  Jtm(3,3,:)---Jxm(4,3,:)-------|----

% Recall that Jxm(i,3,imode+1) and Jtm(i,3,imode+1) are on the radial axis
% Also, Jrm(i,4,imode+1) is half a grid-cell ABOVE r=0 in the prototype

for imode = 0:nmode

% Variable for flipping points below r=0
% NOTE: I THINK SMILEI HAS THIS THE WRONG WAY AROUND:
% PASSING R=0 FOR EVEN M, JX SHOULD NOT CHANGE SIGN, BUT JR AND JT SHOULD
% IF CONSTANT RADIAL JR, THIS POINTS +JY AT THETA=0, -JY AT THETA=PI
% IF CONSTANT RADIAL JX, THIS POINTS +JX AT THETA=0, +JX AT THETA=PI
% CURRENTLY SMILEI IS SETUP TO DO THE OPPOSITE
sign = -1.0;
for i=0:imode+1-1
    sign = sign * -1.0;
end

% Fold Jxm ghost-cells back into simulation window
for i=1:nx+4
    for j=1:2
        Jxm(i,3+j,imode+1) = Jxm(i,3+j,imode+1) + sign * Jxm(i,3-j,imode+1);
        Jxm(i,3-j,imode+1) = sign * Jxm(i,3+j,imode+1);
    end
    if imode > 0
        Jxm(i,3,imode+1) = 0.0;
    else
        % Force dJl/dr = 0 at r=0.
        Jxm(i,3,imode+1) = (4.0*Jxm(i,4,imode+1) - Jxm(i,5,imode+1))/3.0;
   end
end

% Jrm and Jtm boundaries
for i=1:nx+4
    % Fold Jtm
    for j=1:2
        Jtm(i,3+j,imode+1) = Jtm(i,3+j,imode+1) - sign * Jtm(i,3-j,imode+1);
        Jtm(i,3-j,imode+1) = -sign * Jtm(i,3+j,imode+1);
    end

    % Fold Jrm
    for j = 0:1
        Jrm(i,4+j,imode+1) = Jrm(i,4+j,imode+1) - sign * Jrm(i,3-j,imode+1);
        Jrm(i,3-j,imode+1) = -sign * Jrm(i,4+j,imode+1);
    end
    
    % Calculate Jtm and Jrm on r=0 boundary
    if imode == 1
        % Continuity eq. for m=1: div(J) = -drho/dt, but rho=0 on r=0 for
        % m=1, as does Jx. This eq. links Jr and Jt. Also force J_perp to 
        % have no theta dependence: dJ_perp/dtheta = 0. Solve these eq. to
        % get Jtm on r=0
        Jtm(i,3,imode+1) = -Icpx / 8.0 * (9.0 * Jrm(i,4,imode+1) - Jrm(i,5,imode+1));
        % Jrm(3.5) = iJtm(3) from above eq. 
        Jrm(i,3,imode+1) = 2.0*Icpx*Jtm(i,3,imode+1) - Jrm(i,4,imode+1);
    else
        % Jtm=0 for m=0 and m>1, as this implies multiple values at one r
        Jtm(i,3,imode+1) = 0.0;
        % Force Jrm(3.5) = 0, to make Jr=0 on r=0 axis for m=0, m>1
        Jrm(i,3,imode+1) = -Jrm(i,4,imode+1);
    end
end
end
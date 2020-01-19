function [A_1, A_2, B_2] = integration_constants(R, c1, c2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Takes the geometry of the sphere and the concentrations in each material
% to return the integration constants from the boundary conditions:
% 1) zero displacement at r=0
% 2) continuity of radial stress at r=R
% 3) continuity of displacement at r=R
% 4) zero traction at r=1.

% Inputs
% ------
% R  : Nondimensional radius of the core (R_1/R_2)
% c1 : Lithium concentration in the core
% c2 : Lithium concentration in the shell

% Outputs
% -------
% A_1 : First integration constant for the core
% A_2 : First integration constant for the shell
% B_2 : Second integration constant for the shell
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global gamma_1
global gamma_2

%% For simple silicon core, graphite core

[lam1, lam2, G1, G2] = calc_mechanical_constants(c1, c2);

omega = (3.*lam1+2.*G1)*(3.*lam2+2.*G2) + 4*G2*((3.*lam2+2.*G2)*(1-R^3) + (3.*lam1+2.*G1)*R^3); % Eq. 45
A_1 = (1./omega)*((3.*lam1+2.*G1)*((3.*lam2+2.*G2) + 4.*G2*R^3)*gamma_1*c1 + 4.*G2*(1.-R^3)*(3.*lam2+2.*G2)*gamma_2*c2); % Eq. 41
A_2 = (1./omega)*((3.*lam2+2.*G2)*(4.*G2*(1.-R^3) + (3.*lam1+2.*G1))*gamma_2*c2 + 4.*G2*R^3*(3.*lam1+2.*G1)*gamma_1*c1); % Eq. 43
B_2 = (1./omega)*((3.*lam1+2.*G1)*(3.*lam2+2.*G2)*(gamma_1*c1-gamma_2*c2)*R^3); % Eq. 44
end
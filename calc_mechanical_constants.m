function [lam1, lam2, G1, G2] = calc_mechanical_constants(c1, c2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates the lithium-concentration-dependent Lame parameters for both 
% the core and shell materials given by Eq. 49. Here, we calculate the
% nondimensional parameters nondimensionalised by G1 at c1 = 0 given by 
% dim_G_1.
%
% Inputs
% ------
% c1 : Uniform lithium concentration in the core between zero and one, zero
%      is completely unlithiated and one is completely lithiated
% c2 : Uniform lithium concentration in the shell between zero and one, 
%      zero is completely unlithiated and one is completely lithiated
%
% Outputs
% -------
% lam1  : First Lame parameter of the core material
% lam2  : First Lame parameter of the shell material
% G1    : Shear modulus of the core material
% G2    : Shear modulus of the shell material
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global E0_1 E0_2 etaE_1 etaE_2 x_max_1 x_max_2 dim_G_1 nu_1 nu_2

% Calculate the Young's moduli
E_1 = E0_1*(1.+etaE_1*x_max_1*c1); % In Pa as E0_1 given in Pa
E_2 = E0_2*(1.+etaE_2*x_max_2*c2); % In Pa as E0_2 given in Pa

% Calculate the lame parameters and nondimensionalise
lam1 = (((E_1*nu_1)/((1.+nu_1)*(1.-2.*nu_1))))/dim_G_1;
lam2 = (((E_2*nu_2)/((1.+nu_2)*(1.-2.*nu_2))))/dim_G_1;
G1 = (E_1/(2.*(1.+nu_1)))/dim_G_1;
G2 = (E_2/(2.*(1.+nu_2)))/dim_G_1;

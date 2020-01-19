function [lam1, lam2, G1, G2] = calc_mechanical_constants(c1, c2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates the lithiation-dependent Lame parameters for both 
% the core and shell materials given by Eq. 49. Here, we calculate the
% nondimensional parameters nondimensionalised by G1 at c1 = 0 given by 
% dim_G_1.
%
% Inputs
% ------
% c1 : Lithium concentration in the core 
% c2 : Lithium concentration in the shell
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

% Calculate the lame parameters (Eq. 49) and nondimensionalise (Eq. 24)
lam1 = (((E0_1*(1.+etaE_1*x_max_1*c1)*nu_1)/((1.+nu_1)*(1.-2.*nu_1))))/dim_G_1;
lam2 = (((E0_2*(1.+etaE_2*x_max_2*c2)*nu_2)/((1.+nu_2)*(1.-2.*nu_2))))/dim_G_1;
G1 = ((E0_1*(1.+etaE_1*x_max_1*c1))/(2.*(1.+nu_1)))/dim_G_1;
G2 = ((E0_2*(1.+etaE_2*x_max_2*c2))/(2.*(1.+nu_2)))/dim_G_1;

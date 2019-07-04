function constants

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Used to define all the parameters used for the model, including
% mechanical parameters of the core material and the shell material.
% Nondimensional parameters are also calculated. The parameters are then
% made global parameters so the other function can use them.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global fund_charge R_g T
global x_max_1 Jc_SF_1 Vm_1 dim_lam_1 dim_G_1 dim_yield_stress_1 nu_1 E0_1 etaE_1 eta_1 c_max_1 gamma_1 yield_stress_1 Sd_1
global x_max_2 Jc_SF_2 Vm_2 dim_lam_2 dim_G_2 dim_yield_stress_2 nu_2 E0_2 etaE_2 eta_2 c_max_2 gamma_2 yield_stress_2 Sd_2
global c_ratio

%% General Constants
fund_charge = 96483.06883;
R_g = 8.314;
T = 298;

%% Silicon Constants
%%{
x_max_1 = 3.75; % max value of x in Li_xA where A is the core material
Jc_SF_1 = 3.8; % max volume expansion
Vm_1 = 1.2052E-5; % molar volume (for zero lithiation state)
nu_1 = 0.29; % Poisson's Ratio 
E0_1 = 96.E9; % Young's modulus for non-lithiated core material (Pa)
etaE_1 = -0.15278; % Rate which core material's Young's modulus changes per state of charge

eta_1 = (Jc_SF_1 - 1)/(3.*x_max_1); % Coefficient of Compositional Expansion
c_max_1 = x_max_1/Vm_1; % Maximum concentration in mol m^{-3}
%dim_lam_1 = ((E0_1*nu_1)/((1.+nu_1)*(1.-2.*nu_1))); % Dimensional Lame parameter at c=0 in Pa
dim_G_1 = (E0_1/(2.*(1.+nu_1))); % Dimensional Shear Modulus at c=0 in Pa
gamma_1 = 1.; % Non-dimensional concentration pre-factor
%yield_stress_1 = dim_yield_stress_1/(dim_G_1*eta_1*x_max_1); % Non-dimensional yield stress
Sd_1 = (eta_1^2*Vm_1^2*c_max_1*dim_G_1)/(R_g*T); % Stress-assisted diffusion parameter
%}
%% Graphite Constants
%%{
x_max_2 = 1./6.; % max value of x in Li_xB where B is shell material
Jc_SF_2 = 1.1; % max volume expansion
Vm_2 = 8.69E-6; % Molar volume of shell material (for zero lithiation state)
nu_2 = 0.32; % Poisson's Ratio (at c=0)
E0_2 = 32.E9; % Young's modulus for non-lithiated shell material (Pa)
etaE_2 = 14.4375; % Rate which Young's modulus changes per state of charge

eta_2 = (Jc_SF_2 - 1)/(3.*x_max_2); % Coefficient of Compositional Expansion
c_max_2 = x_max_2/Vm_2; % Maximum concentration in mol m^{-3}
%dim_lam_2 = ((E0_2*nu_2)/((1.+nu_2)*(1.-2.*nu_2))); % Dimensional Lame parameter at c=0 in Pa
%dim_G_2 = (E0_2/(2.*(1.+nu_2))); % Dimensional Shear Modulus at c=0 in Pa
gamma_2 = (eta_2*Vm_2*c_max_2)/(eta_1*Vm_1*c_max_1); % Non-dimensional concentration pre-factor
%yield_stress_2 = dim_yield_stress_2/(dim_G_1*eta_1*x_max_1); % Non-dimensional yield stress
Sd_2 = (eta_2*eta_1*Vm_2*Vm_1*c_max_1*dim_G_1)/(R_g*T); % Stress-assisted diffusion parameter
%}

%% Multi-material constants
c_ratio = c_max_2/c_max_1;
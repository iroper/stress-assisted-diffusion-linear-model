function [c1, c2, potential] = SAD_core_shell_conc_calc(R, c0, mu_SF_1, mu_SF_2, options, varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Takes the size of the core, the stress-free chemical potentials of each 
% material, the total amount of lithium in the nano-particle and the 
% options for solving and returns the concentrations in each material and
% the potential throughout the anode.
%
% Inputs
% ------
% R           : Nondimensional radius of the core (R_1/R_2)
% c0          : Fraction of the maximum amount of lithium that is in the 
%               anode particle
% mu_SF_1     : Interpolation function of the non-dimensionalised, 
%               stress-free chemical potential for the core material 
%               against SOC
% mu_SF_2     : Interpolation function of the non-dimensionalised, 
%               stress-free chemical potential for the shell material 
%               against SOC
% options     : Solving options
% plot_option : Option to plot the function that is to be solved
%
% Outputs
% -------
% c1 : Lithium concentration in the core material
% c2 : Lithium concentration in the shell material
% potential : Chemical potential throughout the nano-particle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global c_max_1 gamma_1 Sd_1
global c_max_2 gamma_2 Sd_2
global c_ratio

% If no plot_option provided, assume 0
switch nargin
    case 5
        plot_option = 0;
    case 6
        plot_option = varargin{1};
end

% Calculate bounds to ensure c1 and c2 are between 0 and 1
lb = max(0,(((1.-c0)*c_max_1*R^3)/(c_max_2*(R^3-1))) + c0); % Eq. 48
ub = min(1,c0 - ((c0*c_max_1*R^3)/(c_max_2*(R^3-1)))); % Eq. 48

% Solve for the shell concentration
c2 = solve(@(c2) core_shell_c2_calc(c2,R,mu_SF_1, mu_SF_2, c0), lb, ub, options, plot_option);
% Calculate core concentration given a value of c_0
c1 = c0 + c_ratio*(1-R^(-3))*(c2 - c0); % Eq. 47

% Find the mechanical parameters at these concentrations
[lam1, lam2, G1, G2] = calc_mechanical_constants(c1, c2);

% Calculate potential
if c2>0 && c2<1 % If we are not on the asymptote for the shell OCV
    potential = mu_SF_2(c2) - Sd_2*3.*(3.*lam2 + 2.*G2)*(A_2(c1,c2,R) - gamma_2*c2); % Eq. 35
elseif c1>0 && c1<1 % If we are on the asymptote for shell OCV but not for the core
    potential = mu_SF_1(c1) - Sd_1*3.*(3.*lam1 + 2.*G1)*(A_1(c1,c2,R) - gamma_1*c1); % Eq. 35
else % If both are on asymptote, we cannot calculate the potential
    potential = NaN;
end
end

    function solve_func_eval = core_shell_c2_calc(c2, R, mu_SF_1, mu_SF_2, c0)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Function which evaluates the difference between the chemical
        % potentials of the core and the shell.
        %
        % Inputs
        % ------
        % c2      : Lithium concentration in the shell that is the variable
        %           to be solved for
        % R       : Nondimensional radius of the core (R_1/R_2)
        % mu_SF_1 : Interpolation function of the non-dimensionalised, 
        %           stress-free chemical potential for the core material 
        %           against SOC
        % mu_SF_2 : Interpolation function of the non-dimensionalised, 
        %           stress-free chemical potential for the shell material 
        %           against SOC
        % c0      : Fraction of the maximum amount of lithium that is in
        %           the anode particle
        %
        % Outputs
        % -------
        % solve_func_eval : Residual of equating chemical potentials of 
        %                   the core material and the shell material for a 
        %                   particular value of c2
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        global Sd_1 Sd_2 gamma_1 gamma_2 c_ratio
        c1 = c0 + c_ratio*(1-R^(-3))*(c2 - c0); % compute temporary value of c_Si in function evaluation (Eq. 47)
        % Find the mechanical parameters at the trial concentrations
        [lam1, lam2, G1, G2] = calc_mechanical_constants(c1, c2);
        solve_func_eval = mu_SF_1(c1) - mu_SF_2(c2) - Sd_1*3.*(3.*lam1 + 2.*G1)*(A_1(c1,c2,R) - gamma_1*c1) + Sd_2*3.*(3.*lam2 + 2.*G2)*(A_2(c1,c2,R) - gamma_2*c2); % Eq. 35
    end

    function out = A_1(c1, c2, R)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Function to calculate the core integration constant A_1
        %
        % Inputs
        % ------
        % c1 : Lithium concentration in the core
        % c2 : Lithium concentration in the shell
        % R  : Nondimensional outer radius of the core (R_1/R_2)
        %
        % Outputs
        % -------
        % out : A_1 value
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        global gamma_1 gamma_2
        % Find the mechanical parameters
        [lam1, lam2, G1, G2] = calc_mechanical_constants(c1, c2);
        % Define denominator
        omega = (3.*lam1+2.*G1)*(3.*lam2+2.*G2) + 4*G2*((3.*lam2+2.*G2)*(1-R^3) + (3.*lam1+2.*G1)*R^3); % Eq. 45
        % Calculate A_1 as function of concentrations
        out = (1./omega)*((3.*lam1+2.*G1)*((3.*lam2+2.*G2) + 4.*G2*R^3)*gamma_1*c1 + 4.*G2*(1.-R^3)*(3.*lam2+2.*G2)*gamma_2*c2); % Eq. 41
    end
    
    function out = A_2(c1, c2, R)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Function to calculate shell integration constant A_2
        %
        % Inputs
        % ------
        % c1 : Lithium concentration in the core
        % c2 : Lithium concentration in the shell
        % R  : Nondimensional outer radius of the core (R_1/R_2)
        %
        % Outputs
        % -------
        % out : A_2 value
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        global gamma_1 gamma_2
        % Find the mechanical parameters
        [lam1, lam2, G1, G2] = calc_mechanical_constants(c1, c2);
        % Define denominator
        omega = (3.*lam1+2.*G1)*(3.*lam2+2.*G2) + 4*G2*((3.*lam2+2.*G2)*(1-R^3) + (3.*lam1+2.*G1)*R^3); % Eq. 45
        % Calculate A_2 as function of concentrations
        out = (1./omega)*((3.*lam2+2.*G2)*(4.*G2*(1.-R^3) + (3.*lam1+2.*G1))*gamma_2*c2 + 4.*G2*R^3*(3.*lam1+2.*G1)*gamma_1*c1); % Eq. 43
    end
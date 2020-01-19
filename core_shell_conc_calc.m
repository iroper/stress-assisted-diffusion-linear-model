function [c1, c2, potential] = core_shell_conc_calc(R, c0, mu_SF_1, mu_SF_2, options, varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Takes the volume of the core, the stress-free chemical potentials of each 
% material, the total amount of lithium in the nano-particle and the 
% options for solving and returns the concentrations in each material and
% the potential throughout the anode when stress-assisted diffusion is not
% included in the model.
%
% Inputs
% ------
% R       : Nondimensional radius of the core (R_1/R_2)
% c0      : State of charge
% mu_SF_1 : Interpolation function of the non-dimensionalised, stress-free
%           chemical potential for the core material against SOC
% mu_SF_2 : Interpolation function of the non-dimensionalised, stress-free
%           chemical potential for the shell material against SOC
% options : Solving options
%
% Outputs
% -------
% c1 : Lithium concentration in the core
% c2 : Lithium concentration in the shell
% potential : Chemical potential throughout the particle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global c_max_1
global c_max_2
global c_ratio

% If no plot_option provided, assume 0
switch nargin
    case 5
        plot_option = 0;
    case 6
        plot_option = varargin{1};
end

% Calculate bounds to ensure c_1 and c_2 are between 0 and 1
lb = max(0,(((1.-c0)*c_max_1*R^3)/(c_max_2*(R^3-1))) + c0); % Eq. 48
ub = min(1,c0 - ((c0*c_max_1*R^3)/(c_max_2*(R^3-1)))); % Eq. 48

% Solve for the lithium concentration in the shell
c2 = solve(@(c2) core_shell_c2_calc(c2,R,mu_SF_1, mu_SF_2, c0), lb, ub, options, plot_option);
% Calculate core concentration given a value of c_0
c1 = c0 + c_ratio*(1-R^(-3))*(c2 - c0); % Eq. 47

% Calculate potential
if c2>0 && c2<1 % If we are not on the asymptote for mu_SF_2
    potential = mu_SF_2(c2);
elseif c1>0 && c1<1 % If we are on the asymptote for mu_SF_2 but not for mu_SF_1
    potential = mu_SF_1(c1);
else % If both are on asymptote, we cannot calculate the potential
    disp('Cannot determine potential as on asymptote')
    potential = NaN;
end
end

    function solve_func_eval = core_shell_c2_calc(c2, R, mu_SF_1, mu_SF_2, c0)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Function which evaluates the difference between the two 
        % stress-free chemical potentials of the core and the shell.
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
        % c0      : State of charge
        %
        % Outputs
        % -------
        % solve_func_eval : Residual of equating chemical potentials of 
        %                   the core material and the shell material for a 
        %                   particular value of c2
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        global c_ratio
        c1 = c0 + c_ratio*(1-R^(-3))*(c2 - c0); % compute temporary value of c_1 in function evaluation (Eq. 47)
        solve_func_eval = mu_SF_1(c1) - mu_SF_2(c2);
    end
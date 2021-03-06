function [c1, c2, potential] = SAD_core_shell_conc_calc_special(R, c0, mu_SF_1, mu_SF_2, options, varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Takes the size of the core, the stress-free chemical potentials of each 
% material, the total amount of lithium in the nano-particle and the 
% options for solving and returns the concentrations in each material and
% the potential throughout the anode. This special function evaluates the
% function to be solved for a set of values within the constraints of the
% solution and uses the changes in sign to determine how many roots there
% are. The solver then uses the bisection method close to each change of
% sign to determine each root. If there are multiple roots, the solver
% chooses the lowest one.
%
% Inputs
% ------
% R       : Nondimensional radius of the core (R_1/R_2)
% c0      : Fraction of the maximum amount of lithium that is in the anode
%           particle
% mu_SF_1 : Interpolation function of the non-dimensionalised, stress-free
%           chemical potential for the core material against SOC
% mu_SF_2 : Interpolation function of the non-dimensionalised, stress-free
%           chemical potential for the shell material against SOC
% options : Solving options
%
% Outputs
% -------
% c1 : Lithium concentration in the core material
% c2 : Lithium concentration in the shell material
% potential : Chemical potential throughout the particle
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
c2 = special_solve(@(c2) core_shell_c2_calc(c2,R,mu_SF_1, mu_SF_2, c0), lb, ub, options, plot_option);
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
        % Function which evaluates the difference between the two chemical
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

    function [output, error] = special_solve(solve_fun, lb, ub, options, plot_option)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Function for finding the solution to equations involving the
        % stress-free chemical potentials of the two anode materials. These
        % are tricky as they are numerical functions but need to tend to 
        % infinity at zero and one. Therefore a special solver is required.
        % This solver is different in that it can deal with multiple
        % solutions (in most cases) by evaluating the function to be solved
        % at a range of values between the lower and upper bounds and uses
        % this to determine the number of roots. The bisection method is
        % then used around each change of sign to determine each root.
        %
        % Inputs
        % ------
        % solve_fun   : Function to be solved
        % lb          : Lower bound of the solution
        % ub          : Upper bound of the solution
        % options     : Options for the solver
        % plot_option : Whether the function is to be plotted to check the
        %               solution (0 does not plot and 1 plots)
        %
        % Outputs
        % -------
        % output : Solution that has been found
        % error  : The residual of solve_fun(output)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Produce vector of rough values
        N = 101;
        x_range = linspace(lb,ub,N); % Range of possible solutions
        y_range = zeros(1,N); % Initialise the values of solve_fun
        for c = 1:N
            y_range(c) = solve_fun(x_range(c)); % Evaluate solve_fun at all points in x_range
        end
        if plot_option == 1 % Plot solve_fun
            figure(101)
            plot(x_range,y_range) % Plot
        end
        
        % Find number of solutions
        sol_bounds = find_sol_bounds(x_range,y_range); % Produce matrix of the x values either side of a root
        if size(sol_bounds,1) > 1 % Multiple solutions
            mult_sols = zeros(2,size(sol_bounds,1)); % Intialise the vector of solutions and their errors
            for i = 1:size(sol_bounds,1) % For each change of sign
                % solve for each root using the bisection method
                [mult_sols(1,i),mult_sols(2,i)] = fzero(solve_fun,sol_bounds(i,:));
            end
            % We take the lowest value of the solutions for consistency.
            % The justification around this is if we were charging slowly,
            % it is reasonable to assume that the concentration would only
            % jump to a much higher value once the previous 'branch'
            % disappears.
            output = mult_sols(1,1);
            error = mult_sols(2,1);
        elseif size(sol_bounds,1) > 0 % One solution
            [output,error] = fzero(solve_fun,sol_bounds);
        elseif solve_fun(lb) > 0 && solve_fun(ub) > 0 % All values positive
            % Assuming there is only one intercept with 0, this means all
            % the values of c2 in [lb, ub] yield positive values of the
            % solve_function.
            % However, the difference in the stress-free potentials tends 
            % to infinity at c2 = 1 and to -infinity at c2 = 0.
            % Therefore, if all the function values for c2 in [0,1] are 
            % positive, the only way to get the potentials to equilibrate
            % is if the concentration is at c2=1 and we choose a larger 
            % value of the stress-free potential on the asymptote.
            output = 1.0;
            error = 0;
        elseif solve_fun(lb) < 0 && solve_fun(ub) < 0 % All values negative
            % Same argument as above but for c2 = 0.
            output = 0.0;
            error = 0;
        else % Sometimes, the function at ub or lb is NaN due to the interpolation function not extending all the way to 0 and 1
            [output, error] = fminbnd(@(x) (solve_fun(x)).^2, lb, ub, options);
            
        end
        if error > 10E-4 % If error too big, send to CHECK_ERROR function
            if size(sol_bounds,1)>0
                output = check_error(lb, ub, output, error, solve_fun, options);
                error = 0;
            else
                output = check_error(lb, ub, output, error, @(x) (solve_fun(x)).^2, options);
                error = 0;
            end
        end
        
        if plot_option == 1 % Print solution if the user wishes
            disp('root = ')
            disp(output)
            pause % Allows user to see approximately where the solution is
        end
    end
    
    function sol_bounds = find_sol_bounds(x_range,y_range)
        
        % Find the sign of each function evaluation
        signs = sign(y_range);
        % Find which elements of the vector this changes
        changes_of_sign = find(diff(signs));
        % Trim edges of the vector if there are NaNs at the ends of the
        % signs vector. This occurs when lb != 0 or if ub != 1
        if isnan(signs(1)) == 1
            changes_of_sign = changes_of_sign(2:end);
        end
        if isnan(signs(end)) == 1
            changes_of_sign = changes_of_sign(1:end-1);
        end
        
        sol_bounds = zeros(length(changes_of_sign),2);
        for i = 1:length(changes_of_sign)
            sol_bounds(i,:) = [x_range(changes_of_sign(i)), x_range(changes_of_sign(i)+1)];
        end
    end

    function output = check_error(lb, ub, output, error, optim_fun, options)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Function for finding out why the solution found using SOLVE is
        % above the specified tolerance. As the OCV curves tend to
        % -infinity at 0 and inifity at 1, if the lowest value of the least
        % squares problem is close to the boundary, we assume that the
        % asymptote would take it across zero.
        %
        % Inputs
        % ------
        % lb        : Lower bound of the solution
        % ub        : Upper bound of the solution
        % output    : Solution found using SOLVE
        % error     : Residual of solution found using SOLVE
        % optim_fun : Function to be optimised over
        % options   : Options for the optimiser
        %
        % Outputs
        % -------
        % output : Solution that has been found after finding out why the
        %          residual is high.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if abs(ub-output) < 1E-5
            % If the least-squares solution is very close to the lower
            % bound, we assume the asymptote would make it be the value
            % on the boundary.
            output = ub;
            error = 0; % As we assume this solution is exact, the error is now zero
        elseif abs(lb-output) < 1E-5
            % If the least-squares solution is very close to the upper
            % bound, we assume the asymptote would make it be the value
            % on the boundary.
            output = lb;
            error = 0; % As we assume this solution is exact, the error is now zero
        else % There is some local minimum that the global optimisation has found, therefore, run a global minimiser
            iter = 0;
            while error > 1E-4 && iter < 5
                iter = iter + 1;
                disp('Found local minimum, run global optimisation')
                problem = createOptimProblem('fmincon','objective', optim_fun,'x0',(ub+lb)/2.,'lb',lb,'ub',ub,'options',options);
                ms = MultiStart;
                [output,error] = run(ms,problem,10);
            end
        end
        if error > 10E-4 % If the global minimiser doesn't find a solution with an error under the tolerance, assume no solution
            output = 'There is no solution';
        end
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
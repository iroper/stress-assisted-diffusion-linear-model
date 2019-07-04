function [c1, c2, potential] = core_shell_conc_calc(R, c0, mu_SF_1, mu_SF_2, options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Takes the size of the core, the stress-free chemical potentials of each 
% material, the total amount of lithium in the nano-particle and the 
% options for solving and returns the concentrations in each material and
% the potential throughout the anode when stress-assisted diffusion is not
% included in the model.
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

global c_max_1
global c_max_2
global c_ratio

plot_option = 0; % Plots objective function for debugging

% Calculate bounds to ensure c_Si and c_C are between 0 and 1
lb = max(0,(((1.-c0)*c_max_1*R^3)/(c_max_2*(R^3-1))) + c0);
ub = min(1,c0 - ((c0*c_max_1*R^3)/(c_max_2*(R^3-1))));

% Solve for the lithium concentration in the shell
c2 = solve(@(c2) core_shell_c2_calc(c2,R,mu_SF_1, mu_SF_2, c0), lb, ub, options, plot_option);
% Calculate lithium concentration in the core
c1 = c0 + c_ratio*(1-R^(-3))*(c2 - c0);

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
        % c0      : Fraction of the maximum amount of lithium that is in
        %           the anode particle
        %
        % Outputs
        % -------
        % solve_func_eval : Residual of equating chemical potentials of 
        %                   the core material and the shell material for a 
        %                   particular value of c2
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        global c_ratio
        c1 = c0 + c_ratio*(1-R^(-3))*(c2 - c0); % compute temporary value of c_Si in function evaluation
        solve_func_eval = mu_SF_1(c1) - mu_SF_2(c2);
    end

    function [output, error] = solve(solve_fun, lb, ub, options, plot_option)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Function for finding the solution to equations involving the
        % stress-free chemical potentials of the two anode materials. These
        % are tricky as they are numerical functions but need to tend to 
        % infinity at zero and one. Therefore a special solver is required.
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
        if plot_option == 1
            figure(101)
            % Define x_range and initialise the y_range
            N = 1000;
            x_range = linspace(lb,ub,N); % Range of possible solutions
            y_range = zeros(1,N); % solve_fun for all points in this range
            for c = 1:N
                y_range(c) = solve_fun(x_range(c)); % Evaluate solve_fun at al points in x_range
            end
            plot(x_range,y_range) % Plot
        end
        if solve_fun(lb)*solve_fun(ub) < 0
            % solve_fun has different signs at lb and ub, therefore 
            % solution should be in the interior between these points.
            [output,error] = fzero(solve_fun,[lb,ub]); % solve using the bisection method
            if error > 10E-4 % If the error isn't below a tolerance, send the square of it to the CHECK_ERROR function
                if plot_option == 1
                    disp('error too large')
                end
                output = check_error(lb, ub, output, error, @(x) (solve_fun(x)).^2, options);
                error = 0;
            end
        else % If solve_fun has the same sign, bisection won't work so solve using least squares
            optim_fun = @(x) (solve_fun(x)).^2; % Change solve function to a function to be optimised over
            [output, error] = fminbnd(optim_fun, lb, ub, options);
            if error > 10E-4 % If error too big, send to CHECK_ERROR function
                output = check_error(lb, ub, output, error, optim_fun, options);
                error = 0;
            end
        end
        if plot_option == 1 % Print solution if the user wishes
            disp('c_C = ')
            disp(output)
            pause % Allows user to see approximately where the solution is
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
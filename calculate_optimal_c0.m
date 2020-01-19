function c0 = calculate_optimal_c0(R, mu_SF_1, mu_SF_2, Vmax, options, c0_init, varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to calculate the state of charge which results in an expansion 
% of Vmax, given the size of the core R and the maximum expansion Vmax.
%
% Inputs
% ------
% R       : Nondimensional radius of the core (R_1/R_2)
% mu_SF_1 : Interpolation function of the nondimensional stress-free
%           chemical potential of the core material
% mu_SF_2 : Interpolation function of the nondimensional stress-free
%           chemical potential of the shell material
% Vmax    : The maximum volumetric expansion permitted
% options : Solver options
% c0_init : Initial guess for the state of charge
%
% Outputs
% -------
% c0 : The state of charge which induce an expansion of Vmax in the hybrid
%      particle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    switch nargin
    case 6
        plot_option = 0;
    case 7
        plot_option = varargin{1};
    end
    
    optim_func = @(c0) Vmax_obj_func(c0, R, mu_SF_1, mu_SF_2, Vmax, options);
    [c0, ~] = special_solve_Vmax(optim_func, c0_init, 0.0, 1.0, options, plot_option);
    
end

function res = Vmax_obj_func(c0, R, mu_SF_1, mu_SF_2, Vmax, options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to calculate the residual between the expansion of the
% multi-material particle with an SOC, c0, and a maximum volumetric 
% expansion, Vmax.
%
% Inputs
% ------
% c0      : State of charge
% R       : Nondimensional radius of the core (R_1/R_2)
% mu_SF_1 : Interpolation function of the nondimensional stress-free
%           chemical potential of the core material
% mu_SF_2 : Interpolation function of the nondimensional stress-free
%           chemical potential of the shell material
% Vmax    : The maximum volumetric expansion permitted
% options : Solver options
%
% Outputs
% -------
% res : The difference between the induced expansion and Vmax
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global Vm_1 eta_1 c_max_1

% SAD_core_shell_conc_calc will not work for c_0 not in [0.0,1.0], so if
% the solver tries a c_0 value outside of this range, bring it back to the
% limits.
if c0 > 1.0
    c0 = 1.0;
elseif c0<0.0
    c0 = 0.0;
end

% Calculate the lithium concentration in the core and the shell for this c0
[c1, c2, ~] = SAD_core_shell_conc_calc(R, c0, mu_SF_1, mu_SF_2, options,0);
% Calculating the integration constants of the shell
[~, A_2, B_2] = integration_constants(R, c1, c2);
% Calculating the dimensional displacement at the outer radius
u = eta_1*Vm_1*c_max_1*(A_2 + B_2); % Eq. 59
% Calculate the residual between the relative expanded volume and the
% maximum volume constraint
res = (1+u)^3 - Vmax; % Eq. 59

end

function [output, error] = special_solve_Vmax(solve_fun, c0_init, lb, ub, options, plot_option)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Function for finding the solution to a nonlinear equation in one
    % dimension. This function accounts for the asymptotes in the
    % stress-free chemical potentials at c=0 and c=1. The function has been
    % shown to be multivalued for certain V1 values. It uses the fzero
    % function with an initial guess of the previously calculated c0
    % solution because we want to find the largest solution and this will
    % be the largest because we increase R.
    % Inputs
    % ------
    % solve_fun   : Function to be solved
    % c0_init     : Initial guess for fsolve
    % lb          : Lower bound of the solution
    % ub          : Upper bound of the solution
    % options     : Options for the solver
    % plot_option : Whether the function is to be plotted to check the
    %               solution (0 does not plot and 1 plots)
    %
    % Outputs
    % -------
    % output : Solution that has been found
    % error : The residual of solve_fun(output)
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
        ylim([-0.5,0.5])
        xlim([0,1])
    end
    % Use fsolve with an initial guess
    [output,error] = fzero(solve_fun,c0_init,options);
    % If this solution is invalid, use bisection
    if output < lb || output > ub || abs(error) > 1E-4
        if solve_fun(lb)*solve_fun(ub) < 0
            % solve_fun is different signs at either end therefore solution should be in the interior between these points
            disp('Attempt to find root using fsolve resulted in solution outside constraints, resorting to bisection method')
            % Solve using the bisection method
            [output,error] = fzero(solve_fun,[lb,ub]);
            if abs(error) > 1E-4
                disp('Error still too large, minimising square of function.')
                % Start solve function again using squared function
                [output, error] = special_solve_VV0max(@(x) (solve_fun(x)).^2, c0_init, lb, ub, options, plot_option);
            end
        else
            % Start solve function again using squared function
            disp('Function passes through zero, even number of times. Minimise square')
            [output, error] = special_solve_VV0max(@(x) (solve_fun(x)).^2, c0_init, lb, ub, options, plot_option);
        end
        if error > 1E-4 % If the error isn't below a tolerance, send the square of it to the CHECK_ERROR function
            disp('Error STILL too large. Either on boundary due to asymptotes or global optimiser needed.')
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
        % lb : Lower bound of the solution
        % ub : Upper bound of the solution
        % output : Solution found using SOLVE
        % error : Residual of solution found using SOLVE
        % optim_fun : Function to be optimised over
        % options : Options for the optimiser
        %
        % Outputs
        % -------
        % output : Solution that has been found after finding out why the
        % residual is high.
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


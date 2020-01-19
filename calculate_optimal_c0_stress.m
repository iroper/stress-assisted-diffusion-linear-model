function c0 = calculate_optimal_c0_stress(R, mu_SF_1, mu_SF_2, sigma_max, options, c0_init, varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to calculate the state of charge which results in an effective 
% stress at r=R of sigma_max, given the size of the core R and the maximum
% effective stress sigma_max.
%
% Inputs
% ------
% R         : Nondimensional radius of the core (R_1/R_2)
% mu_SF_1   : Interpolation function of the nondimensional stress-free
%             chemical potential of the core material
% mu_SF_2   : Interpolation function of the nondimensional stress-free
%             chemical potential of the shell material
% sigma_max : The maximum effective stress at r=R permitted
% options   : Solver options
% c0_init   : Initial guess for the optimal state of charge
% varargin  : If plot_option is specified, this determines whether the
%             objective function is plotted
%
% Outputs
% -------
% c0 : The state of charge which induce an effective stress of sigma_max in
%      the hybrid particle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If no plot_option provided, assume 0
switch nargin
    case 6
        plot_option = 0;
    case 7
        plot_option = varargin{1};
end

optim_func = @(c0) sigma_max_obj_func(c0, R, mu_SF_1, mu_SF_2, sigma_max, options);
if plot_option == 1 % Plot the optim_func
    plot_wide(0.0,1.0,optim_func)
end
[c0, ~] = special_solve_sigma_eff_max(optim_func, c0_init, 0.0, 1.0, options, plot_option);

end

function res = sigma_max_obj_func(c0, R, mu_SF_1, mu_SF_2, sigma_max, options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates residual between the von Mises stress at r=R and the maximum
% permitted stress
%
% Inputs
% ------
% c0           : State of charge
% R2           : Dimensional radius of the entire nano-particle
% R            : Nondimensional radius of the core (R_1/R_2)
% mu_SF_1      : Interpolation function of the nondimensional stress-free
%                chemical potential of the core material
% mu_SF_2      : Interpolation function of the nondimensional stress-free
%                chemical potential of the shell material
% sigma_max    : The maximum volumetric expansion permitted
% options      : Solver options
%
% Outputs
% -------
% res : The difference between the induced effective stress at R and 
%       sigma_max
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global Vm_1 eta_1 c_max_1 dim_G_1

% SAD_core_shell_conc_calc will not work for c_0 not in [0.0,1.0]
if c0 > 1.0
    c0 = 1.0;
elseif c0 < 0.0
    c0 = 0.0;
end

% Find c1 and c2 using a special function due to it being multivalued for
% some extreme values of R and c0.
[c1, c2, ~] = SAD_core_shell_conc_calc_special(R, c0, mu_SF_1, mu_SF_2, options, 0);
% Calculate B_2
[~, ~, B_2] = integration_constants(R, c1, c2);
% Calculate G2 for these concentrations
[~, ~, ~, G2] = calc_mechanical_constants(c1, c2);
% Calculate von Mises stress
sigma_eff = (6.*G2*eta_1*Vm_1*c_max_1*dim_G_1*1E-9*abs(B_2))/(R^3); %GPa (Eq. 56)
% Calculate residual
res = sigma_eff - sigma_max;
end

function [output, error] = special_solve_sigma_eff_max(solve_fun, c0_init, lb, ub, options, plot_option)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Function for finding the solution to a nonlinear equation in one
    % dimension. This function accounts for the asymptotes in the
    % stress-free chemical potentials at c=0 and c=1 but also uses a 
    % hard-coded bisection method if the error is too large as this 
    % function has large jumps at points which the fzero bisection method
    % struggles to find.
    %
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
    if solve_fun(lb)*solve_fun(ub) < 0 % solve_fun at each end different signs therefore solution should be in the interior between these points
        % solve using Newton method starting at c0_init
        [output,error] = fzero(solve_fun,c0_init);
        if abs(error) > 10E-4
            % If the error isn't below a tolerance, use bisection instead
            [output,error] = fzero(solve_fun,[lb,ub]);
        end
        if abs(error) > 10E-4
            % If the error still isn't below a tolerance, use hard-coded bisection
            % should return solution to within two floating point numbers
            if plot_option == 1
                [output, ~] = hard_bisection(solve_fun, lb, ub);
                plot_narrow(output,solve_fun)
            else
                [output, ~] = hard_bisection(solve_fun, lb, ub);
            end
            error = 0;
        end
    else % If solve_fun has the same sign, bisection won't work so solve using least squares
        disp('Function same sign at bounds, minimise least squares')
        optim_fun = @(x) (solve_fun(x)).^2; % Change solve function to a function to be optimised over
        if plot_option == 1
            [output, error] = fminbnd(optim_fun, lb, ub, options);
            plot_squared(lb,ub,solve_fun)
        else
            [output, error] = fminbnd(optim_fun, lb, ub, options);
        end
        if error > 10E-4 % If error too big, send to CHECK_ERROR function
            output = check_error(lb, ub, output, error, optim_fun, options);
            error = 0;
        end
    end
    if plot_option == 1
        disp('root = ')
        disp(output)
        pause 
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
            while error > 10E-3 && iter < 5
                iter = iter + 1;
                disp('Found local minimum, run global optimisation')
                problem = createOptimProblem('fmincon','objective', optim_fun,'x0',(ub+lb)/2.,'lb',lb,'ub',ub,'options',options);
                ms = MultiStart;
                [output,error] = run(ms,problem,10);
            end
        end
        if error > 10E-3 % If the global minimiser doesn't find a solution with an error under the tolerance, assume no solution
            output = 'There is no solution';
        end
end

function plot_wide(lb,ub,solve_fun)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plots the function we are trying to solve for all possible values.
    % Includes a zoomed in plot as well to find the root.
    % 
    % Inputs
    % ------
    % lb        : Lower bound for the concentration values to be plotted 
    %             for
    % ub        : Upper bound for the concentration values to be plotted 
    %             for
    % solve_fun : Function to be plotted
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Define x_range and initialise the y_range
    N = 1001;
    x_range = linspace(lb,ub,N); % Range of possible solutions
    y_range = zeros(1,N); % solve_fun for all points in this range
    for c = 1:N
        y_range(c) = solve_fun(x_range(c)); % Evaluate solve_fun at all points in x_range
    end
    
    figure(102)
    subplot(1,2,1)
    plot(x_range,y_range) % Plot whole function
    xlim([lb,ub])
    subplot(1,2,2)
    plot(x_range,y_range) % Plot zoomed in version
    ylim([-0.5,0.5])
    xlim([lb,ub])
end

function plot_narrow(output, solve_fun)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plots the function we are trying to solve for a narrow range of
    % values around the already calculated solution. Includes a zoomed in
    % plot as well.
    % 
    % Inputs
    % ------
    % output    : Solution already calculated
    % solve_fun : Function to be plotted
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Define x_range and initialise the y_range
    N = 1000;
    x_range = linspace(output-0.005,output+0.005,N); % Range of possible solutions around solution
    y_range = zeros(1,N); % solve_fun for all points in this range
    for c = 1:N
        y_range(c) = solve_fun(x_range(c)); % Evaluate solve_fun at all points in x_range
    end
    figure(1012) % figure number reserved for this purpose
    subplot(1,2,1)
    plot(x_range,y_range) % Plot
    subplot(1,2,2)
    plot(x_range,y_range) % Plot
    ylim([-0.5,0.5])
end

function plot_squared(lb,ub,solve_fun)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plots the square of the function we are trying to solve for all
    % possible values of the concentration. Includes a zoomed-in plot.
    % 
    % Inputs
    % ------
    % lb        : Lower bound for the concentration values to be plotted 
    %             for
    % ub        : Upper bound for the concentration values to be plotted 
    %             for
    % solve_fun : Function to be plotted
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Define x_range and initialise the y_range
    N = 1000;
    x_range = linspace(lb, ub, N); % Range of possible solutions
    y_range = zeros(1,N); % solve_fun for all points in this range
    for c = 1:N
        y_range(c) = (solve_fun(x_range(c))).^2; % Evaluate square of the solve_fun at all points in x_range
    end
    figure(1012) % figure number reserved for this purpose
    subplot(1,2,1)
    plot(x_range,y_range) % Plot
    subplot(1,2,2)
    plot(x_range,y_range) % Plot
    ylim([-0.5,0.5])
end

function [x, k] = hard_bisection(func, lb, ub)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculates the root to the function to solve using the bisection
    % method. This should calculate the solution to within two floating
    % point numbers.
    %
    % Inputs
    % ------
    % func : Function to find the root of
    % lb   : Lower bound
    % ub   : Upper bound
    % 
    % Outputs
    % -------
    % x : Solution
    % k : Iterations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    k = 0;
    while abs(ub-lb) > eps*abs(ub) % While the interval is wider than specified
        x = (lb + ub)/2; % test the middle of the interval
        if sign(func(x)) == sign(func(ub)) 
            % If the function at x is the same sign as the function at the upper bound of the interval
            ub = x; % Change upper bound to the midpoint and look in the left half of the interval
        else % Otherwise, look in the right half of the interval
            lb = x;
        end
        k = k + 1; % Increase iteration count
    end
end

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
    % if plot_option is 1, plot the objective function.
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
elseif solve_fun(lb) > 0 && solve_fun(ub) > 0
    % Assuming there is only one intercept with 0, this means all
    % the values of c2 in [lb, ub] yield positive values of the
    % solve_function.
    % However, the difference in the stress-free potentials tends
    % to infinity at c2 = 1 and to -infinity at c2 = 0.
    % Therefore, if all the function values for c2 in [0,1] are
    % positive, the only way to get the potentials to equilibrate
    % is if the concentration is at c2 = 1 and we choose a larger
    % value of the stress-free potential on the asymptote.
    output = 1.0;
    error = 0;
elseif solve_fun(lb) < 0 && solve_fun(ub) < 0
    % Same argument as above but for c2 = 0.
    output = 0.0;
    error = 0;
else % Sometimes, the function at ub or lb is NaN due to the interpolation function not extending all the way to 0 and 1
    optim_fun = @(x) (solve_fun(x)).^2; % Change solve function to a function to be optimised over
    [output, error] = fminbnd(optim_fun, lb, ub, options);
    if error > 10E-4 % If error too big, send to CHECK_ERROR function
        output = check_error(lb, ub, output, error, optim_fun, options);
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
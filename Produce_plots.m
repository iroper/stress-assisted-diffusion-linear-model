%% Driver for producing plots

% Clear previous variables
clear

% Run constants.m to load the constants for silicon, graphite and the
% general constants and make them global
constants;
global fund_charge R_g T
global Vm_1 dim_G_1 eta_1 c_max_1 gamma_1
global gamma_2
global c_ratio

% Set outer radius, options for solvers and linewidth for plots
R2 = 1; % outer radius (used for dimensionalising expansion)
options = optimset('MaxFunEvals', 50000000000000, 'Tolfun', 10E-20, 'MaxIter', 5000000000, 'TolX', 1E-25);
set(0, 'DefaultLineLineWidth', 1.5); % Linewidth increase

% Give core material and shell material
core_material = 'Silicon';
shell_material = 'Graphite';

% Give names of OCV data .mat files
core_OCV_data = 'Silicon_points_interp.mat';
shell_OCV_data = 'Graphite_points_interp.mat';

% Import OCV data and create nondimensional interpolation functions
[OCV_1, OCV_2, mu_SF_1, mu_SF_2] = load_OCV_curves(core_OCV_data, shell_OCV_data);

%% Chemical potentials of individual materials (Figure 2)
%{
figure(2)
% Convert the OCVs to chemical potentials (Eq. 18) and nondimensionalise by
% R_gT (Eq. 20).
plot(OCV_1(:,1), OCV_1(:,2)*((-fund_charge)/(R_g*T)), 'Linewidth', 1.5)
hold on
plot(OCV_2(:,1), OCV_2(:,2)*((-fund_charge)/(R_g*T)), 'Linewidth', 1.5)
ylabel('Stress-free chemical potential, $\mu^{\rm{SF}}_a$','FontSize', 11, 'Interpreter', 'latex')
xlabel('Concentration, $c_a$','FontSize', 11, 'Interpreter', 'latex')
legend({core_material, shell_material}, 'Location', 'East','FontSize', 11, 'Interpreter', 'latex')
set(gca,'FontSize',12)
%}

%% Concentrations, trace and potentials (Figures 3-5)
%{
% Choose c_0 values and organise in vector
c0_N = 99;
c0_min = 0.01;
c0_max = 0.99;
c0_vector = linspace(c0_min, c0_max, c0_N);

% Choose volume of cores
V1_vector = [0.05, 0.1, 0.25, 0.5];
R_vector = V1_vector.^(1./3.);
V1_N = length(V1_vector);

% Initialise concentration and potential matrices
c1_matrix = zeros(V1_N,c0_N);
c2_matrix = zeros(V1_N,c0_N);
pot_matrix = zeros(V1_N,c0_N);
trace_1_matrix = zeros(V1_N, c0_N);
trace_2_matrix = zeros(V1_N, c0_N);

for i = 1:V1_N
    for j = 1:c0_N
        [c1_matrix(i,j),c2_matrix(i,j), pot_matrix(i,j)] = SAD_core_shell_conc_calc(R_vector(i), c0_vector(j), mu_SF_1, mu_SF_2, options, 0);
        [A_1, A_2, B_2] = integration_constants(R_vector(i), c1_matrix(i,j), c2_matrix(i,j));
        [lam1, lam2, G1, G2] = calc_mechanical_constants(c1_matrix(i,j), c2_matrix(i,j));
        trace_1_matrix(i,j) = 3.*(3.*lam1+2.*G1)*(A_1 - gamma_1*c1_matrix(i,j)); % Eq. 36 + 2*Eq. 37
        trace_2_matrix(i,j) = 3.*(3.*lam2+2.*G2)*(A_2 - gamma_2*c2_matrix(i,j)); % Eq. 36 + 2*Eq. 37
    end
end

color_vec = ['b', 'k', 'r', 'g', 'c', 'm', 'b', 'k', 'r', 'g'];
% Plot concentrations (Figure 3)
figure(3)
plot(c0_vector, c0_vector, color_vec(1))
legend_labels = cell(1,2*V1_N+1);
legend_labels(1) = {'Single Material'};
for i = 1:V1_N
    hold on
    plot(c0_vector, c1_matrix(i,:), strcat('-',color_vec(i+1)))
    legend_labels(2*i) = {['$\psi = ' num2str(V1_vector(i)) '$ $(c_1)$']};
    hold on
    plot(c0_vector, c2_matrix(i,:), strcat('--',color_vec(i+1)))
    legend_labels(2*i+1) = {['$\psi = ' num2str(V1_vector(i)) '$ $(c_2)$']};
end
ylabel('Concentration in each material, $c_1$ \& $c_2$','FontSize', 16, 'Interpreter', 'latex')
xlabel('State of charge, $c_0$','FontSize', 16, 'Interpreter', 'latex')
legend(legend_labels, 'Location', 'Southeast','FontSize', 10, 'Interpreter', 'latex')
set(gca,'FontSize',12)

% Plot traces (Figure 4)
figure(4)
plot(c0_vector, zeros(1,length(c0_vector)), color_vec(1))
legend_labels = cell(1,2*V1_N+1);
legend_labels(1) = {'Single Material'};
for i = 1:V1_N
    hold on
    plot(c0_vector, trace_1_matrix(i,:), strcat('-',color_vec(i+1)))
    legend_labels(2*i) = {['$\psi = ' num2str(V1_vector(i)) '$ $(\Omega_1)$']};
    hold on
    plot(c0_vector, trace_2_matrix(i,:), strcat('--',color_vec(i+1)))
    legend_labels(2*i+1) = {['$\psi = ' num2str(V1_vector(i)) '$ $(\Omega_2)$']};
end
ylabel('Trace of Cauchy Stress Tensor, tr$(\boldmath{\sigma})$','FontSize', 16, 'Interpreter', 'latex')
xlabel('State of charge, $c_0$','FontSize', 16, 'Interpreter', 'latex')
legend(legend_labels, 'Location', 'Southwest','FontSize', 10, 'Interpreter', 'latex')
set(gca,'FontSize',12)

% Plot potentials (Figure 5)
figure(5)
plot(OCV_2(:,1), OCV_2(:,2)*((-fund_charge)/(R_g*T)), color_vec(1))
legend_labels = cell(1,V1_N+2);
legend_labels(1) = {shell_material};
for i = 1:V1_N
    hold on
    plot(c0_vector, pot_matrix(i,:), color_vec(i+1))
    legend_labels(i+1) = {['$\psi = ' num2str(V1_vector(i)) '$']};
end
hold on
plot(OCV_1(:,1), OCV_1(:,2)*((-fund_charge)/(R_g*T)), color_vec(V1_N+2))
legend_labels(end) = {core_material};
ylabel('Chemical potential, $\mu$','FontSize', 16, 'Interpreter', 'latex')
xlabel('State of charge, $c_0$','FontSize', 16, 'Interpreter', 'latex')
legend(legend_labels, 'Location', 'Northwest','FontSize', 14, 'Interpreter', 'latex')
set(gca,'FontSize',12)
%}

%% Concentration and potentials without stress-assisted diffusion (Figures 6 & 7)
%{
% Choose c_0 values and organise in vector
c0_N = 99;
c0_min = 0.01;
c0_max = 0.99;
c0_vector = linspace(c0_min, c0_max, c0_N);

% Choose volumes of the core
V1_vector = [0.05, 0.1, 0.25, 0.5];
R_vector = V1_vector.^(1./3.);
V1_N = length(V1_vector);

% Initialise concentration and potential matrices
c1_matrix = zeros(V1_N,c0_N);
c2_matrix = zeros(V1_N,c0_N);
pot_matrix = zeros(V1_N,c0_N);

for i = 1:V1_N
    for j = 1:c0_N
        [c1_matrix(i,j),c2_matrix(i,j), pot_matrix(i,j)] = core_shell_conc_calc(R_vector(i), c0_vector(j), mu_SF_1, mu_SF_2, options);
    end
end

color_vec = ['b', 'k', 'r', 'g', 'c', 'm', 'b', 'k', 'r'];

% Plot concentrations (Figure 6)
figure(6)
plot(c0_vector, c0_vector, color_vec(1))
legend_labels = cell(1,2*V1_N+1);
legend_labels(1) = {'Single Material'};
for i = 1:V1_N
    hold on
    plot(c0_vector, c1_matrix(i,:), strcat('-',color_vec(i+1)))
    legend_labels(2*i) = {['$\psi = ' num2str(V1_vector(i)) '$ $(c_1)$']};
    hold on
    plot(c0_vector, c2_matrix(i,:), strcat('--',color_vec(i+1)))
    legend_labels(2*i+1) = {['$\psi = ' num2str(V1_vector(i)) '$ $(c_2)$']};
end
ylabel('Concentration in each material, $c_1$ and $c_2$','FontSize', 16, 'Interpreter', 'latex')
xlabel('State of charge, $c_0$','FontSize', 16, 'Interpreter', 'latex')
legend(legend_labels, 'Location', 'Northwest','FontSize', 10, 'Interpreter', 'latex')
set(gca,'FontSize',12)

% Plot potentials (Figure 7)
figure(7)
plot(OCV_2(:,1), OCV_2(:,2)*((-fund_charge)/(R_g*T)), color_vec(1))
legend_labels = cell(1,V1_N+2);
legend_labels(1) = {shell_material};
for i = 1:V1_N
    hold on
    plot(c0_vector, pot_matrix(i,:), color_vec(i+1))
    legend_labels(i+1) = {['$\psi = ' num2str(V1_vector(i)) '$']};
end
hold on
plot(OCV_1(:,1), OCV_1(:,2)*((-fund_charge)/(R_g*T)), color_vec(V1_N+2))
legend_labels(end) = {core_material};
ylabel('Chemical potential, $\mu$','FontSize', 16, 'Interpreter', 'latex')
xlabel('State of charge, $c_0$','FontSize', 16, 'Interpreter', 'latex')
legend(legend_labels, 'Location', 'Southeast','FontSize', 14, 'Interpreter', 'latex')
set(gca,'FontSize',12)
ylim([-20,0]) % For aesthetics
%}

%% V, Q/V and sigma_eff (Figures 8-10)
%{
% Choose c_0 values and organise in vector
c0_N = 9;
c0_min = 0.1;
c0_max = 0.9;
c0_vector = linspace(c0_min, c0_max, c0_N);

% Choose volumes of the core
V1_vector = linspace(0.01,0.99,99);
R_vector = V1_vector.^(1./3.);
V1_N = length(V1_vector);

% Define dimensionalisation factor of displacement
disp_dim_factor = R2*eta_1*Vm_1*c_max_1; % Eq. 20
stress_dim_factor = eta_1*Vm_1*c_max_1*dim_G_1*1E-9; %GPa % Eq. 20

% Initialise V, Q/V and sigma_eff matrices
V_vec = zeros(c0_N+1,V1_N);
QV_vec = zeros(c0_N+1,V1_N);
sigma_eff_vec = zeros(c0_N+1,V1_N);

for j = 1:V1_N
    for i = 1:c0_N
        [c1, c2, potential] = SAD_core_shell_conc_calc(R_vector(j), c0_vector(i), mu_SF_1, mu_SF_2, options, 0);
        [~, ~, ~, G2] = calc_mechanical_constants(c1, c2);
        [A_1, A_2, B_2] = integration_constants(R_vector(j), c1, c2);
        [u, ~, ~, ~, ~] = calculate_displacement(A_1, A_2, B_2, R_vector(j), 1);
        V_vec(i,j) = (1.+disp_dim_factor*u)^3; % Eq. 52
        sigma_eff_vec(i,j) = (4.*G2*stress_dim_factor*abs(B_2))/((R_vector(j))^3); % Eq. 56
        QV_vec(i,j) = (c1*(R_vector(j))^3 + c2*c_ratio*(1-(R_vector(j))^3))/((1.+disp_dim_factor*u)^3); % Eq. 57
    end
    % Compute for c_0 = 1
    [A_1, A_2, B_2] = integration_constants(R_vector(j), 1.0, 1.0);
    [u, u1, u2, r1, r2] = calculate_displacement(A_1, A_2, B_2, R_vector(j), 1);
    [~, ~, ~, G2] = calc_mechanical_constants(1.0, 1.0);
    V_vec(c0_N+1,j) = (1.+u*disp_dim_factor)^3; % Eq. 52
    sigma_eff_vec(c0_N+1,j) = (4.*G2*stress_dim_factor*abs(B_2))/((R_vector(j))^3); % Eq. 56
    QV_vec(c0_N+1,j) = (1.0*(R_vector(j))^3 + 1.0*c_ratio*(1-(R_vector(j))^3))/((1.+u*disp_dim_factor)^3); % Eq. 57
end

colour_vec = {'blue', 'black', 'red', 'green', 'cyan', 'magenta', 'blue', 'black', 'red', 'green'}';
line_style_vec = {'-', '-', '-', '-', '-', '-', '-.', '-.', '-.', '-.'}';

% Plot V (Figure 8)
figure(8)
legend_labels = cell(1,c0_N+1);
for i = 1:c0_N
    h(i) = plot(V1_vector, V_vec(i,:), 'Linewidth', 1.5);
    hold on
    legend_labels(i) = {['$c_0 = ' num2str(c0_vector(i)) '$']};
end
h(i+1) = plot(V1_vector, V_vec(end,:), 'Linewidth', 1.5);
legend_labels(end) = {'$c_0 = 1.0$'};
hold off
xlabel('Core volume fraction, $\psi$','FontSize', 16, 'Interpreter', 'latex')
ylabel('Relative expanded volume, $V$','FontSize', 16, 'Interpreter', 'latex')
set(h, {'Color'}, colour_vec);
set(h, {'LineStyle'}, line_style_vec);
legend(h,legend_labels, 'Location', 'Northwest','FontSize', 12, 'Interpreter', 'latex')
set(gca,'FontSize',12)

% Plot Q/V (Figure 9)
figure(9)
legend_labels = cell(1,c0_N+1);
for i = 1:c0_N
    h(i) = plot(V1_vector, QV_vec(i,:), 'Linewidth', 1.5);
    hold on
    legend_labels(i) = {['$c_0 = ' num2str(c0_vector(i)) '$']};
end
h(i+1) = plot(V1_vector, QV_vec(end,:), 'Linewidth', 1.5);
legend_labels(end) = {'$c_0 = 1.0$'};
hold off
xlabel('Core volume fraction, $\psi$','FontSize', 16, 'Interpreter', 'latex')
ylabel({'Relative amount of lithium per'; 'expanded volume, $Q/V$'},'FontSize', 12, 'Interpreter', 'latex')
set(h, {'Color'}, colour_vec);
set(h, {'LineStyle'}, line_style_vec);
set(gca,'FontSize',12)
% Use columnlegend and ylim for aesthetics
columnlegendQV(4, legend_labels, 'Location', 'northwest','FontSize', 11.5, 'Interpreter', 'latex', 'boxon');
ylim([0.0,0.27])
% Use generic legend
%legend(legend_labels, 'Location', 'best','FontSize', 10, 'Interpreter', 'latex')

% Plot sigma_eff (Figure 10)
figure(10)
legend_labels = cell(1,c0_N+1);
for i = 1:c0_N
    h(i) = plot(V1_vector, sigma_eff_vec(i,:), 'Linewidth', 1.5);
    hold on
    legend_labels(i) = {['$c_0 = ' num2str(c0_vector(i)) '$']};
end
h(i+1) = plot(V1_vector, sigma_eff_vec(end,:), 'Linewidth', 1.5);
legend_labels(end) = {'$c_0 = 1.0$'};
hold off
xlabel('Core volume fraction, $\psi$','FontSize', 16, 'Interpreter', 'latex')
ylabel('Maximum induced stress, $\sigma^*_{\rm{eff}}(R_{\rm{1}})$, (GPa)','FontSize', 12, 'Interpreter', 'latex')
set(h, {'Color'}, colour_vec);
set(h, {'LineStyle'}, line_style_vec);
set(gca,'FontSize',12)
% Use columnlegend and ylim for aesthetics of paper
columnlegendsigma(2, legend_labels, 'Location', 'northwest','FontSize', 12, 'Interpreter', 'latex', 'boxon');
ylim([0.0,110])
% Use generic legend
%legend(h, legend_labels, 'Location', 'southeast','FontSize', 10, 'Interpreter', 'latex')
%}

%% Q_max for prescribed V_max (Figure 11)
%{

plot_concs = 1; % Change to one to plot concentrations
plot_option = 0; % Change to one to debug the solve function

% Choose volumes of core
V1_vector_origin = linspace(0.01,0.99,99);
R_vector_origin = V1_vector_origin.^(1./3.);
V1_N = length(V1_vector_origin);

% Choose V_max values to test
V_max_vector = [1.2, 1.4, 1.6, 1.8];
V_max_N = length(V_max_vector);

% Define eta_bar_1
eta_bar_1 = eta_1*Vm_1*c_max_1;

% Intialise Q_max matrix and concentration vectors if they are wanted
Q_matrix = zeros(V_max_N,V1_N+1);
if plot_concs == 1
    c1_matrix = zeros(V_max_N,V1_N+1);
    c2_matrix = zeros(V_max_N,V1_N+1);
    c0_matrix = zeros(V_max_N,V1_N+1);
end

% Initialise different V1_vectors
V1_matrix = zeros(V_max_N, V1_N+1);

% Need c0_init for the solver as there are multiple solutions for some
% values of V1 and thus need to find the maximum c_0 solution, and
% c_0 seems to be monotonically decreasing, thus using the previous 
% solution as a starting point is a good idea.
c0_init = 1.0;

for i = 1:V_max_N % For each expansion constraint
    % Calculate R_hat using Eq. A.2
    [lam1, lam2, G1, G2] = calc_mechanical_constants(1.0, 1.0);
    R_lim1 = ((3.*lam1+2.*G1)*(3.*lam2+2.*G2) + 4.*G2*(3.*lam2+2.*G2))*((V_max_vector(i))^(1./3.) - 1) - eta_bar_1*(3.*lam2+2.*G2)*(4.*G2+3.*lam1+2.*G1)*gamma_2; 
    R_lim2 = eta_bar_1*3.*(3.*lam1+2.*G1)*(lam2+2.*G2)*gamma_1 - eta_bar_1*(3.*lam2+2.*G2)*(4.*G2+3.*lam1+2.*G1)*gamma_2 - 4.*G2*((3.*lam1+2.*G1) - (3.*lam2+2.*G2))*((V_max_vector(i))^(1./3.) - 1);
    R_lim = (R_lim1/R_lim2)^(1./3.);
    % Check if this value is physical using Eq. A.3 
    lower_bnd = (((eta_bar_1*(3.*lam2 + 2.*G2)*gamma_2*(3.*lam1 + 2.*G1 + 4.*G2))/((3.*lam1+2.*G1)*(3.*lam2+2.*G2) + 4.*G2*(3.*lam2+2.*G2))) + 1.)^3;
    if V_max_vector(i)>lower_bnd && V_max_vector(i)<(1.+eta_bar_1*gamma_1)^3
        % If so, add this R_lim to R_vector_origin to avoid a kink in the plot
        R_vector = [R_vector_origin(R_vector_origin<R_lim), R_lim, R_vector_origin(R_vector_origin>R_lim)];
    else % If not, repeat first entry so that the vectors are the right length
        R_vector = [R_vector_origin(1), R_vector_origin];
    end
    V1_matrix(i,:) = R_vector.^3; % Update the V_1 values with new addition
    
    for j = 1:V1_N+1 % For each size of silicon core
        if R_vector(j) <= R_lim
            % If the constraint is not reached at c_0 = 1.0, Q, c_1 and 
            % c_2 are easily calculated 
            Q_matrix(i,j) = 1.0*(R_vector(j))^3 + 1.0*c_ratio*(1-(R_vector(j))^3); % Eq. 51
            if plot_concs == 1
                c1_matrix(i,j) = 1.0;
                c2_matrix(i,j) = 1.0;
                c0_matrix(i,j) = 1.0;
            end
            c0_init = 1.0;
        else
            % If we need to constrain the lithium concentration below 1.0,
            % solve for the value of c_0 which gives V = V_max
            c0_opt = calculate_optimal_c0(R2, R_vector(j), mu_SF_1, mu_SF_2, V_max_vector(i), options, c0_init, plot_option);
            % Calculate the c1 and c2 values at this c0 value
            [c1_opt, c2_opt, ~] = SAD_core_shell_conc_calc(R_vector(j), c0_opt, mu_SF_1, mu_SF_2, options, 0);
            
            if plot_concs == 1
                c1_matrix(i,j) = c1_opt;
                c2_matrix(i,j) = c2_opt;
                c0_matrix(i,j) = c0_opt;
            end
            
            % Calculate the charge inside the microparticle
            Q_matrix(i,j) = c1_opt*(R_vector(j))^3 + c2_opt*c_ratio*(1-(R_vector(j))^3); % Eq. 51
            % Update the initial guess
            c0_init = c0_opt;
        end
    end
end

color_vec = ['b', 'k', 'r', 'g', 'c', 'm', 'b', 'k', 'r', 'g'];

% Plot Q_max against the core volume (Figure 11)
figure(11)
legend_labels = cell(1,V_max_N);
for i = 1:V_max_N
    plot(V1_matrix(i,:), Q_matrix(i,:), color_vec(i))
    hold on
    legend_labels(i) = {['$V_{\rm{max}} = ' num2str(V_max_vector(i)) '$']};
end
ylabel('Maximum amount of lithium, $Q_{\rm{max}}$','FontSize', 16, 'Interpreter', 'latex')
xlabel('Core volume fraction, $\psi$','FontSize', 16, 'Interpreter', 'latex')
legend(legend_labels, 'Location', 'Northeast','FontSize', 12, 'Interpreter', 'latex')
set(gca,'FontSize',12)

if plot_concs == 1
    for i = 1:V_max_N
        figure(111+i)
        plot(V1_matrix(i,:), c1_matrix(i,:), 'r')
        hold on
        plot(V1_matrix(i,:), c2_matrix(i,:), 'b')
        hold on
        plot(V1_matrix(i,:), c0_matrix(i,:), 'm')
        title(['$V_{\rm{max}} = ' num2str(V_max_vector(i)) '$'],'FontSize', 16, 'Interpreter', 'latex')
        ylabel('Concentrations','FontSize', 16, 'Interpreter', 'latex')
        xlabel('Core volume fraction, $\psi$','FontSize', 16, 'Interpreter', 'latex')
        legend({'$c_1$','$c_2$', '$c_0$'}, 'FontSize', 14, 'Interpreter', 'latex')
    end
end
%}

%% Q_max for prescribed maximum stress (Figures 12-13)
%{
% Choose V1 values
V1_vector = linspace(0.01,0.99,99);
R_vector_origin = V1_vector.^(1./3.);
V1_N = length(V1_vector);

% Choose V1 values for the sigma_eff vs c_0 plot
V1_vector2 = linspace(0.1,0.9,9);
R_vector2 = V1_vector2.^(1./3.);
V_N2 = length(V1_vector2);
c0_vector2 = linspace(0.01,0.99,99);
c0_N2 = length(c0_vector2);
sigma_matrix = zeros(V_N2,c0_N2+1);

% Choose sigma_max values to test
sigma_max_vector = [0.5, 1.0, 2.0, 3.0, 4.0, 5.0]; %GPa
sigma_max_N = length(sigma_max_vector);

% Define dimensionalisation factor of displacement
disp_dim_factor = R2*eta_1*Vm_1*c_max_1;
stress_dim_factor = eta_1*Vm_1*c_max_1*dim_G_1*1E-9; %GPa

% Intialise Q_maxvector
Q_matrix = zeros(sigma_max_N,V1_N+1);
V1_matrix = zeros(sigma_max_N, V1_N+1);

plot_concs = 0; % Change to one to plot concentrations
plot_option = 0; % Change to one to debug the solve function

if plot_concs == 1
    c1_matrix = zeros(sigma_max_N,V1_N+1);
    c2_matrix = zeros(sigma_max_N,V1_N+1);
    c0_matrix = zeros(sigma_max_N,V1_N+1);
end

for i = 1:V_N2
    for j = 1:c0_N2
        [c1, c2, ~] = SAD_core_shell_conc_calc(R_vector2(i), c0_vector2(j), mu_SF_1, mu_SF_2, options, 0);
        [~, ~, ~, G2] = calc_mechanical_constants(c1, c2);
        [~, ~, B_2] = integration_constants(R_vector2(i), c1, c2);
        sigma_matrix(i,j) = (6.*G2*eta_1*Vm_1*c_max_1*dim_G_1*1E-9*abs(B_2))/(R_vector2(i)^3); % Eq. 56
    end
    [~, ~, ~, G2] = calc_mechanical_constants(1.0, 1.0);
    [~, ~, B_2] = integration_constants(R_vector2(i), 1.0, 1.0);
    sigma_matrix(i,end) = (6.*G2*eta_1*Vm_1*c_max_1*dim_G_1*1E-9*abs(B_2))/(R_vector2(i)^3); % Eq. 56
end

colour_vec = {'blue', 'black', 'red', 'green', 'cyan', 'magenta', 'blue', 'black', 'red'}';
line_style_vec = {'-', '-', '-', '-', '-', '-', '-.', '-.', '-.'}';

% Plot the maximum induced stress against c_0 (Figure 12)
figure(12)
legend_labels = cell(1,V_N2);
for i = 1:V_N2
    w(i) = plot([c0_vector2, 1.0],sigma_matrix(i,:));
    hold on
    legend_labels(i) = {['$\psi = ' num2str(V1_vector2(i)) '$']};
end
set(w, {'Color'}, colour_vec);
set(w, {'LineStyle'}, line_style_vec);
xlabel('State of charge, $c_0$','FontSize', 16, 'Interpreter', 'latex')
ylabel('Maximum induced stress, $\sigma^*_{\rm{eff}}(R_1)$, (GPa)','FontSize', 12, 'Interpreter', 'latex')
legend(w,legend_labels, 'Location', 'Southeast','FontSize', 12, 'Interpreter', 'latex')
set(gca,'FontSize',12)

% Plot the inset figure in Figure 12
figure(125)
for i = 1:V_N2
    w(i) = plot([c0_vector2,1.0],sigma_matrix(i,:));
    hold on
end
set(w, {'Color'}, colour_vec);
set(w, {'LineStyle'}, line_style_vec);
set(gca,'FontSize',12)
grid on
% These are the limits which work for the example in the paper, they will
% need to be changed for other parameters
ylim([0.0,20.0])
xlim([0.05,0.25])

% Initiate initial guess
c0_init = 1.0;
for i = 1:sigma_max_N % For each expansion constraint
    % Calculate R_hat which give V = V_max for c_0 = 1.0 using Eq. A.5
    [lam1, lam2, G1, G2] = calc_mechanical_constants(1.0, 1.0);
    R_lim1 = (3.*lam1+2.*G1)*(3.*lam2+2.*G2)*(6.*stress_dim_factor*G2*abs(gamma_1-gamma_2) - sigma_max_vector(i)) - 4.*G2*(3.*lam2+2.*G2)*sigma_max_vector(i);
    R_lim2 = 4.*G2*((3.*lam1+2.*G1) - (3.*lam2+2.*G2))*sigma_max_vector(i);
    V1_lim = R_lim1/R_lim2;
    
    % Lower bound for sigma_max values which allow c_0=1 for any V1
    saturated_lb = (6.*stress_dim_factor*G2*(3.*lam1+2.*G1)*(3.*lam2+2.*G2)*abs(gamma_1-gamma_2))/((3.*lam1+2.*G1)*(3.*lam2+2.*G2) + 4.*G2*(3.*lam2+2.*G2));
    % Upper bound for sigma_max values which allow c_0=1 for any V1
    saturated_ub = (6.*stress_dim_factor*G2*(3.*lam2+2.*G2)*abs(gamma_1-gamma_2))/(3.*(lam2+2.*G2));
    
    R_lim = (V1_lim)^(1./3.);
    % Check if this value is physical
    if V1_lim>0 && V1_lim<1 % If so, add this to the R values to avoid a kink in the plot
        R_vector = [R_vector_origin(R_vector_origin<R_lim), R_lim, R_vector_origin(R_vector_origin>R_lim)];
    else % If not, repeat first entry so that the vectors are the right length
        R_vector = [R_vector_origin(1), R_vector_origin];
    end
    V1_matrix(i,:) = R_vector.^3; % Update the V1 values with new addition
    
    for j = 1:V1_N+1 % For each size of silicon core
        if V1_lim>0 && V1_lim<1 && R_vector(j) <= R_lim
            % If the constraint is not reached at c_0 = 1.0, Q, c1 and c2
            % are easily calculated 
            Q_matrix(i,j) = 1.0*(R_vector(j))^3 + 1.0*c_ratio*(1-(R_vector(j))^3);
            c0_init = 1.0;
            if plot_concs == 1
                c1_matrix(i,j) = 1.0;
                c2_matrix(i,j) = 1.0;
                c0_matrix(i,j) = 1.0;
            end
        else
            % Solve for the value of c_0 which gives V = V_max
            c0_opt = calculate_optimal_c0_stress(R_vector(j), mu_SF_1, mu_SF_2, sigma_max_vector(i), options, c0_init, plot_option);

            % Calculate the c1 and c2 values at this c_0 value
            [c1_opt, c2_opt, ~] = SAD_core_shell_conc_calc_special(R_vector(j), c0_opt, mu_SF_1, mu_SF_2, options,0);
            if plot_concs == 1
                c1_matrix(i,j) = c1_opt;
                c2_matrix(i,j) = c2_opt;
                c0_matrix(i,j) = c0_opt;
            end

            % Calculate the amount of Li inside the nano-particle
            Q_matrix(i,j) = c1_opt*(R_vector(j))^3 + c2_opt*c_ratio*(1-(R_vector(j))^3); % Eq. 51
        end
    end
end

color_vec = ['b', 'k', 'r', 'g', 'c', 'm', 'b', 'k', 'r', 'g'];

% Plot Q_max for a maximum induced stress (Figure 13)
figure(13)
legend_labels = cell(1,sigma_max_N);
for i = 1:sigma_max_N
    plot(V1_matrix(i,:), Q_matrix(i,:), color_vec(i))
    hold on
    legend_labels(i) = {['$\sigma_{\rm{max}} = ' num2str(sigma_max_vector(i)) '$']};
end
ylabel('Maximum amount of lithium, $Q_{\rm{max}}$','FontSize', 16, 'Interpreter', 'latex')
xlabel('Core volume fraction, $\psi$','FontSize', 16, 'Interpreter', 'latex')
set(gca,'FontSize',12)
columnlegendsigma(2, legend_labels, 'Location', 'northwest','FontSize', 12, 'Interpreter', 'latex', 'boxon');
ylim([0.0,0.09])
%legend(legend_labels, 'Location', 'best','FontSize', 14, 'Interpreter', 'latex')

if plot_concs == 1
    for i = 1:sigma_max_N
        figure(13+i)
        plot(V1_matrix(i,:), c1_matrix(i,:), 'r')
        hold on
        plot(V1_matrix(i,:), c2_matrix(i,:), 'b')
        hold on
        plot(V1_matrix(i,:), c0_matrix(i,:), 'm')
        title(['$\sigma_{\rm{max}} = ' num2str(sigma_max_vector(i)) '$'],'FontSize', 16, 'Interpreter', 'latex')
        ylabel('Concentrations','FontSize', 16, 'Interpreter', 'latex')
        xlabel('Core volume fraction, $\psi$','FontSize', 16, 'Interpreter', 'latex')
        legend({'$c_1$','$c_2$', '$c_0$'}, 'FontSize', 14, 'Interpreter', 'latex')
    end
end

% Note that the graphite concentration for large sigma_max has some erratic
% behaviour for large volumes of the silicon core. This is because for
% large R and small c_0, there are cases where there is a repeated root of
% for c2. Due to the numerical solving, it is hard to determine whether 
% these are true roots or not, hence the oscillation between cases when the
% solver determines that the root is true and those which it determines 
% them to not be true roots.
% At these low c2 levels and the small volume of graphite, the variations
% do not change the optimal value of c_0 by very much at all, and thus our
% result is a very good approximation, hence we ignore these fluctuations.

%}
















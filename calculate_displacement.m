function [u, u1, u2, r1, r2] = calculate_displacement(A_1, A_2, B_2, R, r)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Takes the integration constants for core and the shell and calculates
% the displacement within each material using equation 34.

% Inputs
% ------
% A_1  : Integration constant for the core
% A_2  : First integration constant for shell
% B_2  : Second integration constant for shell
% R    : Nondimensional radius of the core (R_1/R_2)
% r    : Radius or radii that the displacement is to be calculated for

% Outputs
% -------
% u  : Vector of the displacement at the radial points within the whole of 
%        the vector r
% u1 : Vector of the displacement at the radial points within Omega_1
% u2 : Vector of the displacement at the radial points within Omega_2
% r1 : Vector of the radial points within the Omega_1
% r2 : Vector of the radial points within the Omega_2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if R > 1 % This is only valid for cores and shells thus R must be less than 1
    disp('Error, R must be greater than 1')
elseif max(r)>1 || min(r) < 0 % Radial points must be within [0, 1]
    disp('Error, value(s) of r to be calculated must be between zero and 1')
elseif max(r)<= R % If all the radial points are within the Omega_1
    r1 = r;
    r2 = [];
    u1 = A_1*r1;
    u2 = [];
    u = u1;
elseif min(r) >= R % If all the radial points are within Omega_2
    r1 = [];
    r2 = r;
    u1 = [];
    u2 = A_2*r2 + (B_2./(r2.^2));
    u = u2;
else % If the radial points cross over both regions
    if isempty(r(abs(r-R)<1E-10)) % If there is no radial point corresponding to R
        disp('Warning: R is not an element of the vector r. This will cause there to be an incorrectly interpolated solution')
        % Warning of the kink that will occur in the plot of the displacement profile
    end
    % Separate the radial points, keeping R within each vector
    r1 = r(r<=R);
    r2 = r(r>=R);
    % Calculate the displacement within each material
    u1 = A_1*r1;
    u2 = A_2*r2 + (B_2./(r2.^2));
    % Concatenate the displacement profiles, only allowing one R point.
    u = [A_1*(r(r<=R)),A_2*(r(r>R)) + (B_2./((r(r>R)).^2))];
end
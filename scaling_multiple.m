clear all; clc;

% ====================================
% Planet and Core Constants
% ====================================
R = 6052e3;      % Input Planet Radius [m]
R_core = 3110e3; % Input Planet core radius [m]

A_planet = pi * R^2;
A_core = pi * R_core^2;
A_annulus = A_planet - A_core;

V_planet = (4/3) * pi * R^3;
V_core   = (4/3) * pi * R_core^3;
V_annulus = V_planet - V_core;

% ====================================
% Load Input Data
% ====================================
data_3D = load('impacts_3D.dat');
num_rows = size(data_3D, 1);
data_2D = zeros(num_rows, 4);

% Copy 1st and 3rd columns unchanged
data_2D(:, 1) = data_3D(:, 1);
data_2D(:, 3) = data_3D(:, 3);

% Set all values in the 4th column to 1.570796
data_2D(:, 4) = 1.570796;

% ====================================
% Loop Over Each Impactor Radius
% ====================================
for i = 1:num_rows
    r_imp = data_3D(i, 2);
    r_sphere = (3^(1/3)) * r_imp;  % Input Isocore Radius for 3D

    % Calculate 3D isocore volume
    d_3D  = R - r_sphere/2;
    x0_3D = (R^2 - r_sphere^2 + d_3D^2) / (2 * d_3D);

    V1 = pi * ( ((x0_3D * r_sphere^2) - ((x0_3D^3)/3) + d_3D*x0_3D^2 - x0_3D*d_3D^2) - ...
                ((d_3D-r_sphere)*r_sphere^2 - (((d_3D-r_sphere)^3)/3) + d_3D*(d_3D-r_sphere)^2 - (d_3D-r_sphere)*d_3D^2) );

    V2 = pi * ( R^3 - ((R^3)/3) - x0_3D*R^2 + ((x0_3D^3)/3) );

    V = V1 + V2;
    V_ratio = V / V_annulus;

    % Solve for 2D r_circle
    objective_func = @(r_c) (calc_area_2D(r_c, R) / A_annulus) - V_ratio;
    r_circle_v1 = sqrt(r_sphere^3 * (R^2 - R_core^2) / (R^3 - R_core^3));
    r_circle_v2 = fzero(objective_func, r_circle_v1);

    % Scaled radius of the impactor (2nd column)
    r_imp_scaled_v2 = r_circle_v2 / (3^(1/3));
    data_2D(i, 2) = r_imp_scaled_v2;
end

% ====================================
% Write Output File (impacts_2D.dat)
% ====================================
fid = fopen('impacts_2D.dat', 'w');
for i = 1:num_rows
    fprintf(fid, '  %8.6f  %11.2f    %8.6f    %8.6f\n', ...
        data_2D(i, 1), data_2D(i, 2), data_2D(i, 3), data_2D(i, 4));
end
fclose(fid);
disp('impacts_2D.dat successfully generated.');

% ================= FUNCTION =================
function A = calc_area_2D(r_c, R)
    d  = R - r_c/2;
    x0 = (R^2 - r_c^2 + d^2) / (2 * d);

    % A1 is the total area under the curve of the small circle (Isocore)
    A1 = (x0 - d)*sqrt(r_c^2 - (x0 - d)^2) + (r_c^2)*asin((x0 - d)/r_c) + (pi * r_c^2)/2;

    % A2 is the total area under the curve of the big circle (Planet)
    A2 = (pi * R^2)/2 - x0*sqrt(R^2 - x0^2) - (R^2)*asin(x0/R);

    A = A1 + A2;
end

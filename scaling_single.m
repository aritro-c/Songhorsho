clear all; clc;

R = 6052e3; % Input Planet Radius [m]
R_core = 3110.e3; % Input Planet core radius [m]
r_imp = 157106.04; % Impactor radius [m]

r = (3^(1/3))*r_imp;  % Input Isocore Radius for 3D


% r_circle = ?


% ====================================
% Planet and Core Areas and Volume:
% ====================================

A_planet = pi*R^2;
A_core = pi*R_core^2;
A_annulus = A_planet - A_core;


V_planet = (4/3) * pi * R^3;     
V_core   = (4/3) * pi * R_core^3; 
V_annulus = V_planet - V_core;



% ===============================
% Calculate 3D isocore volume
% ===============================
r_sphere = r;
d_3D  = R - r_sphere/2;
x0_3D = (R^2 - r_sphere^2 + d_3D^2) / (2 * d_3D);


% V1 is the total volume under the curve of the small circle (Isocore)
V1 = pi * ( (   (x0_3D * r_sphere^2) - ((x0_3D^3)/3) + d_3D*x0_3D^2 - x0_3D*d_3D^2   )  -  ((d_3D-r_sphere)*r_sphere^2 - (((d_3D-r_sphere)^3)/3) + d_3D*(d_3D-r_sphere)^2 -(d_3D-r_sphere)*d_3D^2 ) );

% V2 is the total volume under the curve of the big circle (Planet)
V2 = pi*(  R^3 - ((R^3)/3) - x0_3D*R^2 +  ((x0_3D^3)/3) );


% Total volume, V = V1 + V2
V = V1+V2;
V_old = (4/3)*pi*r_sphere^3;
V_ratio = V/V_annulus;

% ====================================
% Solve for 2D r_circle
% ====================================

% Find at what r_circle value A_ratio == V_ratio
objective_func = @(r_c) (calc_area_2D(r_c, R) / A_annulus) - V_ratio;

% Initial guess from Version-1 (full circle/sphere approximation)
r_circle_v1 = sqrt(r_sphere^3 * (R^2 - R_core^2) / (R^3 - R_core^3));

% Find exact r_circle for Version-2
r_circle_v2 = fzero(objective_func, r_circle_v1);

% Verify that A_ratio matches V_ratio
A_final = calc_area_2D(r_circle_v2, R);
A_ratio = A_final / A_annulus;


r_imp_scaled_v1 = r_circle_v1/(3^(1/3));
r_imp_scaled_v2 = r_circle_v2/(3^(1/3));

r_imp/1000
r_imp_scaled_v1
r_imp_scaled_v2



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



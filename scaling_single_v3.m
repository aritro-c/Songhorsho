clear all; clc;

R = 6052e3; % Input Planet Radius [m]
R_core = 3110.e3; % Input Planet core radius [m]
r_imp = 521951.61; % Impactor radius [m]

r = (3^(1/3))*r_imp;  % Isocore Radius


% r_spherical_annulus = ?


% ====================================
% Planet and Core Volumes:
% ====================================


% 3D
V_planet_3D = (4/3) * pi * R^3; 
V_planet_core_3D = (4/3) * pi * R_core^3; 
V_mantle_3D = V_planet_3D - V_planet_core_3D;

% 2D Spherial Annulus
V_mantle_annulus = ((2*pi)/3)*((R^3)-(R_core^3));


% ====================================
% 3D Isocore Volume:
% ====================================

r_sphere = r;
d_3D  = R - r_sphere/2;
x0_3D = (R^2 - r_sphere^2 + d_3D^2) / (2 * d_3D);


% V1 is the total volume of the isocore inside the curve of the small circle (Isocore)
V1_isocore_3D = pi * ( (   (x0_3D * r_sphere^2) - ((x0_3D^3)/3) + d_3D*x0_3D^2 - x0_3D*d_3D^2   )  -  ((d_3D-r_sphere)*r_sphere^2 - (((d_3D-r_sphere)^3)/3) + d_3D*(d_3D-r_sphere)^2 -(d_3D-r_sphere)*d_3D^2 ) );

% V2 is the total volume of the isocore inside the curve of the big circle (Planet)
V2_isocore_3D = pi*(  R^3 - ((R^3)/3) - x0_3D*R^2 +  ((x0_3D^3)/3) );


% Total volume, V = V1 + V2

V_isocore_3D = V1_isocore_3D+V2_isocore_3D;

V_ratio_3D = V_isocore_3D/V_mantle_3D;


% ====================================
% Solve for r_spherical_annulus
% ====================================

% Find at what r_spherical_annulus value V_ratio_2D == V_rato_3D
objective_func = @(r_ic) (calc_volume_spherical_annulus(r_ic,R)/V_mantle_annulus) - V_ratio_3D;

% Initial guess from Version-1 (full circle/sphere approximation)
r_circle_v1 = sqrt(r_sphere^3 * (R^2 - R_core^2) / (R^3 - R_core^3));

% Find exact r_spherical_annulus for Version-3
r_spherical_annulus = fzero(objective_func, r_circle_v1);

% Verify that V_ratio_2D matches V_ratio_3D
V_2D = calc_volume_spherical_annulus(r_spherical_annulus,R);
V_ratio_2D = V_2D/V_mantle_annulus;


% Find out the scaled radius of the impactor
r_imp_scaled_v1 = r_circle_v1/(3^(1/3));
r_imp_scaled_v3 = r_spherical_annulus/(3^(1/3));

r_imp
r_imp_scaled_v1
r_imp_scaled_v3




% ================= FUNCTION =================
function V = calc_volume_spherical_annulus(r_ic, R)

    d = R-r_ic/2;
    x_0 = (R^2 - r_ic^2 + d^2) / (2 * d);

    %% Calculate the volume inside the isocore curve:

    % First, define the function that represents the "area" of the vertical strip:
    S_x_1 = @(x) sqrt((x.^2 + r_ic.^2 - (x-d).^2).*(r_ic.^2 - (x-d).^2)) + (x.^2).*log((sqrt(x.^2 + r_ic.^2 - (x-d).^2) + sqrt(r_ic.^2 - (x-d).^2))./x);

    % Calculate the definite integral:
    V1 = integral(S_x_1,(d-r_ic),x_0); 


    %% Calculate the volume inside the planet curve:

    % First, define the function that represents the "area" of the vertical strip:
    S_x_2 = @(x) R.*sqrt(R.^2 - x.^2)  +  (x.^2).*log((R + sqrt(R.^2 - x.^2))./x);


    % Calculate the definite integral:
    V2 = integral(S_x_2,(x_0),R);

    %% Total volume of the isocore:
    V = V1 + V2;

end

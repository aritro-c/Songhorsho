clear all; clc;



R = 6000e3;  % Input Planet Radius
R_core = 1000e3; % Input Planet core radius

r = 100e3;  % Input Isocore Radius







% ==============================================================

R_circle = R; 
R_core_circle = R_core;
r_circle = r; 

R_sphere = R;
R_core_sphere = R_core;
r_sphere = r; 

% ============================================================

% 2D
d_circle = R_circle - r_circle/2; % This is the distance of the center of the isocore from the center of the planet (from StagYY module advdif.f90)
x0_circle = (R_circle^2-r_circle^2+d_circle^2)/(2*d_circle); % This is the value of the X-coordinate where Planet and isocore intersect, i.e., the limit for the integration 


% 3D
d_sphere = R_sphere - r_sphere/2; % This is the distance of the center of the isocore from the center of the planet (from StagYY module advdif.f90)
x0_sphere = (R_sphere^2-r_sphere^2+d_sphere^2)/(2*d_sphere); % This is the value of the X-coordinate where Planet and isocore intersect, i.e., the limit for the integration 


% =============================================================

A_planet = pi*R_circle^2; 
A_core = pi*R_core_circle^2;

A_annulus = A_planet - A_core;



V_planet = (4/3)*pi*R_sphere^3;
V_core = (4/3)*pi*R_core_sphere^3;

V_annulus = V_planet - V_core;




% ==============================================================


% A1 is the total area under the curve of the small circle (Isocore)
A1 = ( (x0_circle-d_circle)*sqrt(r_circle^2 - (x0_circle-d_circle)^2) + (r_circle^2)*(asin((x0_circle-d_circle)/(r_circle))) + (pi*r_circle^2)/(2));

% A2 is the total area under the curve of the big circle (Planet)
A2 = (pi*R_circle^2)/2 - x0_circle*sqrt(R_circle^2 - x0_circle^2) - R_circle^2*asin(x0_circle/R_circle);


% Total area, A = A1 + A2

A = A1+A2;

A_old = pi*r_circle^2;


% ===================================================================



% V1 is the total volume under the curve of the small circle (Isocore)

V1 = pi * ( (   (x0_sphere * r_sphere^2) - ((x0_sphere^3)/3) + d_sphere*x0_sphere^2 - x0_sphere*d_sphere^2   )  -  ((d_sphere-r_sphere)*r_sphere^2 - (((d_sphere-r_sphere)^3)/3) + d_sphere*(d_sphere-r_sphere)^2 -(d_sphere-r_sphere)*d_sphere^2 ) );

% V2 is the total volume under the curve of the big circle (Planet)

V2 = pi*(  R_sphere^3 - ((R_sphere^3)/3) - x0_sphere*R_sphere^2 +  ((x0_sphere^3)/3) );


% Total volume, V = V1 + V2

V = V1+V2;

V_old = (4/3)*pi*r_sphere^3;


% ==================================================================


A_ratio = A/A_annulus;
V_ratio = V/V_annulus;


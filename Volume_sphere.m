clear all; clc;


R = 5; % Input Planet Radius
r = 3; % Input Isocore Radius


d = R - r/2; % This is the distance of the center of the isocore from the center of the planet (from StagYY module advdif.f90)
x0 = (R^2-r^2+d^2)/(2*d); % This is the value of the X-coordinate where Planet and isocore intersect, i.e., the limit for the integration 


% V1 is the total volume under the curve of the small circle (Isocore)

V1 = pi * ( (   (x0 * r^2) - ((x0^3)/3) + d*x0^2 - x0*d^2   )  -  ((d-r)*r^2 - (((d-r)^3)/3) + d*(d-r)^2 -(d-r)*d^2 ) );

% A2 is the total area under the curve of the big circle (Planet)

V2 = pi*(  R^3 - ((R^3)/3) - x0*R^2 +  ((x0^3)/3) );


% Total volume, V = V1 + V2

V = V1+V2;

V_old = (4/3)*pi*r^3;

clear all; clc;

R = 5; % Input Planet Radius
r = 3; % Input Isocore Radius




d = R - r/2; % This is the distance of the center of the isocore from the center of the planet (from StagYY module advdif.f90)
x0 = (R^2-r^2+d^2)/(2*d); % This is the value of the X-coordinate where Planet and isocore intersect, i.e., the limit for the integration 

% A1 is the total area under the curve of the small circle (Isocore)
A1 = ( (x0-d)*sqrt(r^2 - (x0-d)^2) + (r^2)*(asin((x0-d)/(r))) + (pi*r^2)/(2));

% A2 is the total area under the curve of the big circle (Planet)
A2 = (pi*R^2)/2 - x0*sqrt(R^2 - x0^2) - R^2*asin(x0/R);


% Total area, A = A1 + A2

A = A1+A2;

A_old = pi*r^2;

clear all; clc;

R_annulus = 5; % Input Planet Radius
r_annulus = 3; % Input Isocore Radius

d = R_annulus - r_annulus/2;  % This is the distance of the center of the isocore from the center of the planet (from StagYY module advdif.f90)
x_0 = ((R_annulus^2)-(r_annulus^2)+d^2)/(2*d); % This is the value of the X-coordinate where Planet and isocore intersect, i.e., the limit for the integration 



%% Calculate the volume inside the isocore curve:

% First, define the function that represents the "area" of the vertical strip:
S_x_1 = @(x) sqrt((x.^2 + r_annulus.^2 - (x-d).^2).*(r_annulus.^2 - (x-d).^2)) + (x.^2).*log((sqrt(x.^2 + r_annulus.^2 - (x-d).^2) + sqrt(r_annulus.^2 - (x-d).^2))./x);

% Calculate the definite integral:
V1 = integral(S_x_1,(d-r_annulus),x_0); 





%% Calculate the volume inside the planet curve:

% First, define the function that represents the "area" of the vertical strip:
S_x_2 = @(x) R_annulus.*sqrt(R_annulus.^2 - x.^2)  +  (x.^2).*log((R_annulus + sqrt(R_annulus.^2 - x.^2))./x);


% Calculate the definite integral:
V2 = integral(S_x_2,(x_0),R_annulus);


%% Total volume of the isocore:
V = V1 + V2;

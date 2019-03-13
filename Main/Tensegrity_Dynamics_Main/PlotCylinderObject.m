function [LatFace, UpFace, DwFace] = PlotCylinderObject(p1,p2,r,nsurfpatches)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 

% Returns points to generate surface plots for the 
% lateral (LatFace), upper (UpFace), and down (DwFace) of a cylinder
% LatFace.x, LatFace.y, LatFace.z: Coordinates of lateral cylinder surface
% Same for UpFace, DwFace
%
% p1: Start point of cylinder
% p2: End point of cylinder
% r: Radius of cylinder
% nsurfpatches: # number of surface patches for the cylinder


eta1 = linspace(0,1,nsurfpatches+1);
eta2 = linspace(0,1,2);
[eta1, eta2] = meshgrid(eta1, eta2);

dir = p2 - p1; % Direction of cylinder

%% Calculating orthogonal vectors to cylinder surface
[maxn,maxi] = max(abs(dir));
    
switch maxi
    case 1
        v1 = [-dir(2)/dir(1); 1; 0]/norm([-dir(2)/dir(1); 1; 0],2);
        v2 = [-dir(3)/dir(1); 0; 1]/norm([-dir(3)/dir(1); 0; 1],2); 
    case 2
        v1 = [1; -dir(1)/dir(2); 0]/norm([1; -dir(1)/dir(2); 0],2);
        v2 = [0; -dir(3)/dir(2); 1]/norm([0; -dir(3)/dir(2); 1],2);         
    case 3
        v1 = [1; 0; -dir(1)/dir(3)]/norm([1; 0; -dir(1)/dir(3)],2);
        v2 = [0; 1; -dir(2)/dir(3)]/norm([0; 1; -dir(2)/dir(3)],2);
end

%% Lateral face
LatFace.x = p1(1) + r.*cos(2.*pi.*eta1).*v1(1) ...
    + r.*sin(2.*pi.*eta1).*v2(1) + dir(1).*eta2;
LatFace.y = p1(2) + r.*cos(2.*pi.*eta1).*v1(2) ...
    + r.*sin(2.*pi.*eta1).*v2(2) + dir(2).*eta2;
LatFace.z = p1(3) + r.*cos(2.*pi.*eta1).*v1(3) ...
    + r.*sin(2.*pi.*eta1).*v2(3) + dir(3).*eta2;

%% Upper face
UpFace.x = p1(1) + eta2.*r.*cos(2.*pi.*eta1).*v1(1) ...
    + eta2.*r.*sin(2.*pi.*eta1).*v2(1);
UpFace.y = p1(2) + eta2.*r.*cos(2.*pi.*eta1).*v1(2) ...
    + eta2.*r.*sin(2.*pi.*eta1).*v2(2);
UpFace.z = p1(3) + eta2.*r.*cos(2.*pi.*eta1).*v1(3) ...
    + eta2.*r.*sin(2.*pi.*eta1).*v2(3);

%% Down face
DwFace.x = p2(1) + eta2.*r.*cos(2.*pi.*eta1).*v1(1) ...
    + eta2.*r.*sin(2.*pi.*eta1).*v2(1);
DwFace.y = p2(2) + eta2.*r.*cos(2.*pi.*eta1).*v1(2) ...
    + eta2.*r.*sin(2.*pi.*eta1).*v2(2);
DwFace.z = p2(3) + eta2.*r.*cos(2.*pi.*eta1).*v1(3) ...
    + eta2.*r.*sin(2.*pi.*eta1).*v2(3);
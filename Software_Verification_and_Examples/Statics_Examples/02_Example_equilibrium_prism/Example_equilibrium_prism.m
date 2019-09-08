%% Prism minimum mass calculation, solve NK=W, NP=D
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 
%
% EXAMPLE:
%    - Manual node matrix specification
%    - Manual connectivity matrix generation
%    - Prestress equilibrium solver
%    - Pinned_nodes = [n1 1 0 1; n2 1 0 0];
%      n1 node fixed in x and z direction, n2 node fixed in x direction 

% % EXAMPLE:
clear all; clc; close all; warning off;
% 
% Specify node positions
N = 1*[0.5 0 0; 0 0.866 0; -0.5 0 0; 0.5 0 1; 0 0.866 1; -0.5 0 1]';

% Specify bar connectivity
Cb_in = [3 5; 1 6; 2 4]; % Bar 1 connects node 3 to 5, etc
C_b = tenseg_ind2C(Cb_in,N);

% Specify string connectivity
Cs_in = [1 2; 2 3; 3 1; 4 5; 5 6; 6 4; 1 4; 2 5; 3 6]; % String one is node 1 to 2
C_s = tenseg_ind2C(Cs_in,N);

% Plot structure before string segmentation
tenseg_plot(N,C_b,C_s);

%% ------------------------ Other Option

prism.Bradius = 0.02*ones(size(C_b,1),1); % Radius of bars [# bars x 1]
prism.Sradius = 0.01*ones(size(C_s,1),1); % Radius of strings [# strings x 1]
prism.Nradius = 0.005*ones(size(N,2),1); % Radius of node spheres [# nodes x 1]

tenseg_plot(N,C_b,C_s,[],[],[],[],prism)

%% Another Prism Configuration
% R = 6; z=20;
% N = [-sqrt(3)/2*R sqrt(3)/2*R 0 1/sqrt(3)*R -1/(2*sqrt(3))*R -1/(2*sqrt(3))*R;
%      -1/2*R -1/2*R 1*R 0 1/2*R -1/2*R;
%      0 0 0 z z z];
% 
% %  Specify bar connectivity
% Cb_in = [3 6; 1 4; 2 5]; % Bar 1 connects node 3 to 5, etc
% C_b = tenseg_ind2C(Cb_in,N);
% 
% Specify string connectivity
% Cs_in = [1 2; 2 3; 3 1; 4 5; 5 6; 6 4; 3 4; 1 5; 2 6]; % String one is node 1 to 2
% C_s = tenseg_ind2C(Cs_in,N);


%% Perform the calculation
W = zeros(size(N));
% W(3,1:3) = 1000*ones(1,3);
W(3,4:6) = -1000000*ones(1,3);

prism.N = N;
prism.C_b = C_b;
prism.C_s = C_s;
prism.W = W;
prism.Pinned_nodes = [4 1 1 1; 5 1 1 1; 6 1 1 1];
prism.bar_material='Aluminum'; % Specify bar material Aluminum, UHMWPE or Steel
prism.string_material='Aluminum'; % Specify bar material Aluminum, UHMWPE or Steel
prism.bar_failure='yielding_n_buckling'; % gives maximum mass subject to yielding and buckling

% Solve for an equilibrium condition
[Force_den,Support_Load,LABEL,MIN_MASS,Mass_bar,Mass_string,BarWidth_in,StringWidth_in,Loop] = tenseg_equilibrium_minimal_mass(prism);

BarWidth=BarWidth_in;
StringWidth=StringWidth_in;
Loop=Loop-1;
tenseg_plot_physical(N,C_b,C_s,[],[],[],[],'Sradius',LABEL,BarWidth,StringWidth);
% title('D Bar Structure')

save Example_equilibrium_prism

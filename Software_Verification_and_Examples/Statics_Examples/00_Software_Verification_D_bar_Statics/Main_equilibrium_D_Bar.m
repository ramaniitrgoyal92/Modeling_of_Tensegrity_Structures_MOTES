%% D-Bar minimum mass calculation, solve NK=W, NP=D
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 

% EXAMPLE:
%    - Manual node matrix specification
%    - Manual connectivity matrix generation
%    - Prestress equilibrium solver
%    - Pinned_nodes = [n1 1 0 1; n2 1 0 0];
%      n1 node fixed in x and z direction, n2 node fixed in x direction
clear all; clc; close all;
% Manually specify node positions (in 3D).
n1 = [-1 0 0]';
n2 = [0 -1 0]';
n3 = [1 0 0]';
n4 = [0 1 0]';
% Put node vectors in node matrix. Node matrix has to be 3xn for n nodes.
N = [n1 n2 n3 n4];
% Manually specify connectivity indices.
C_s_in = [1 3;   % This is indicating that string 1 is the vector from node 1 to node 3,
    2 4];  %    and that bar 2 connects node 2 to node 4.
C_b_in = [1 2;   % Similarly, this is saying bar 1 connects node 1 to node 2,
    2 3;   %    and so on...
    3 4;
    4 1];
% Convert the above matrices into full connectivity matrices.
C_b = tenseg_ind2C(C_b_in,N);
C_s = tenseg_ind2C(C_s_in,N);
% Plot the structure to make sure it looks right
tenseg_plot(N,C_b,C_s);
% title('D Bar Structure')
%%
% segments = 3; % Divide strings into a number of segments
% [N,C_b,C_s,parents] = tenseg_string_segment(N,C_b,C_s,segments);
% tenseg_plot(N,C_b,C_s); % Plot structure to see what string segmentation looks like

%%

% Define the Force Matrix
W = zeros(size(N));
W(1,1) = 10000; W(1,3) = -10000;
% Pinned_nodes = [2 1 1 1; 4 1 1 1]; % Node 2, in x and z direction
Pinned_nodes = []; % Node 2, in x and y direction
% Define a tensegrity data structure, named 'kite'
D_Bar.N = N;
D_Bar.C_b = C_b;
D_Bar.C_s = C_s;
D_Bar.W = W;
D_Bar.Pinned_nodes = Pinned_nodes;
D_Bar.bar_material='Aluminum'; % Specify bar material Aluminum, UHMWPE or Steel
D_Bar.string_material='Aluminum'; % Specify bar material Aluminum, UHMWPE or Steel
D_Bar.bar_failure='yielding_n_buckling'; % gives maximum mass subject to yielding and buckling
% % Solve for an equilibrium condition
% [Force_den,Support_Load] = tenseg_equilibrium_V1(D_Bar);
% % Solve for mass of bars and strings
% [total_mass,bar_mass,string_mass] = Mass_statics(D_Bar,Force_den);

% Solve for an equilibrium condition
[Force_den,Support_Load,LABEL,MIN_MASS,Mass_bar,Mass_string,BarWidth_in,StringWidth_in,Loop] = tenseg_equilibrium_minimal_mass(D_Bar);

BarWidth=BarWidth_in;
StringWidth=StringWidth_in;
Loop=Loop-1;
tenseg_plot_physical(N,C_b,C_s,[],[],[],[],'Sradius',LABEL,BarWidth,StringWidth);
% title('D Bar Structure')

save D_Bar_Statics.mat

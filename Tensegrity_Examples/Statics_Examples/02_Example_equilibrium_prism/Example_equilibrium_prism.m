% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 

clear all; close all; clc

[N,C_b,C_s] = tenseg_prismplate(1,1,1);
tenseg_plot(N,C_b,C_s);

prism.N = N;
prism.C_b = C_b;
prism.C_s = C_s;

W = zeros(size(N));
% W(3,1:3) = -1000*ones(1,3);
% W(3,4:6) = 1*ones(1,3);

prism.W = W;
prism.Pinned_nodes = [4 1 1 1; 5 1 1 1; 6 1 1 1];


% Solve for an equilibrium condition 
[Force_den,Support_Load] = tenseg_equilibrium(prism);

% save Example_equilibrium_prism.mat
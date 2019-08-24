% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 

% EXAMPLE:
%    - Manual node matrix specification
%    - Manual connectivity matrix generation
%    - Prestress equilibrium solver
%    - Pinned_nodes = [n1 1 0 1; n2 1 0 0]; n1 node fixed in x and z direction, n2 node fixed in x direction 

clear all; clc; close all; warning off

% Manually specify node positions (in 3D).
n1 = [-1 1/2 0]';
n2 = [-1/2 -1 0]';
n3 = [1 1/2 0]';
n4 = [-1/2 1 0]';

% Put node vectors in node matrix. Node matrix has to be 3xn for n nodes.
N = [n1 n2 n3 n4];

% Manually specify connectivity indices.
C_b_in = [1 3;   % This is indicating that bar 1 is the vector from node 1 to node 3,
	      2 4];  %    and that bar 2 connects node 2 to node 4.

C_s_in = [1 2;   % Similarly, this is saying string 1 connects node 1 to node 2,
	      2 3;   %    and so on...
		  3 4;
		  4 1]; 

% Convert the above matrices into full connectivity matrices.
C_b = tenseg_ind2C(C_b_in,N);
C_s = tenseg_ind2C(C_s_in,N);

% Plot the structure to make sure it looks right
tenseg_plot(N,C_b,C_s);

% Define the Force Matrix
W = zeros(size(N));
W(:,1) = [1 0 0];
W(:,3) = [-1 0 0];
W(:,2) = [0 1 0];
W(:,4) = [0 -1 0];


Pinned_nodes = [2 1 1 1;3 1 1 1];

% Define a tensegrity data structure, 'kite'
kite.N = N;
kite.C_b = C_b;
kite.C_s = C_s;
kite.W = W;
kite.Pinned_nodes = Pinned_nodes;

% Solve for an equilibrium condition 
[Force_den,Support_Load] = tenseg_equilibrium(kite);

% save Example_equilibrium_kite.mat

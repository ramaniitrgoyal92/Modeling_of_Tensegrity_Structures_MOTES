% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 

% EXAMPLE:
clear all; clc; close all;


% Specify initial node positions for class k structure

%   - Note that you don't have to worry about converting class k nodes into
%     constrained class 1 physical and virtual nodes yourself.

n1 = [0 0 0]';
n2 = [1 0 0]';
n3 = [1 1 0]';
n4 = [0 1 0]';
N_simple = [n1 n2 n3 n4];

% Specify bar and string connectivity 
C_b_in = [4 1;  % Bar 1 connects node 4 to node 1
	      1 2;  % Bar 2 connects node 1 to node 2
		  2 3;  % Bar 3 connects node 2 to node 3
		  3 4]; % Bar 4 connects node 3 to node 4
C_s_in = [2 4;  % String 1 connects node 2 to node 4
	      1 3]; % String 2 connects node 1 to node 3

% Convert above index notation into actual connectivity matrices
C_b = tenseg_ind2C(C_b_in,N_simple);
C_s = tenseg_ind2C(C_s_in,N_simple);

% pinned_nodes
pinned_nodes = [1];

% To illustrate the class k conversion process, print the initial number of
% nodes here
disp(['Initial # of nodes: ' num2str(size(N_simple,2))])
tenseg_plot(N_simple,C_b,C_s);

% Convert specified class k structure into a class 1 structure with constraints
[N_new,C_b_new,C_s_new,P,D,node_constraints] = tenseg_class_k_convert(N_simple,C_b,C_s,pinned_nodes);


% Print the final number of nodes
disp(['Converted Class K # of nodes:' num2str(size(N_new,2))])

% Print the generated node constraints
disp(['Class K Node constraints:'])
for i=1:numel(node_constraints)
	if numel(node_constraints{i})>1
		disp(['   Coincident nodes: ' num2str(node_constraints{i})])
	end
end

%Specify resting string lengths

% Here, we're setting every string rest length to 70% of its given length
S_0_percent = [(1:size(C_s_new,1))', 0.7*ones(size(C_s_new,1),1)]; % percent of initial lengths

% This function converts those specified percentages into rest lengths
s_0 = tenseg_percent2s0(N_new,C_s_new,S_0_percent);



%% Simulation 

% Create data structure of system BEFORE segmentation
classK_test.N = N_new;
classK_test.C_b = C_b_new;
classK_test.C_s = C_s_new;
classK_test.P = P;
classK_test.D = D;
classK_test.s_0 = [1;0.5];
classK_test.tf = 5;
classK_test.dt = .001;
classK_test.video = 1;

%%
% Perform simulation
[History,debug] = tenseg_sim_classkopen(classK_test);

%% Plotting position and velocity of specified node/axes
% tenseg_plot_node(Hist,[1,2],[1,2,3])
% tenseg_plot_velocity(Hist,[1,2],[1,2,3])

%% Create plots with 3D objects 
classK_test.Bradius = [0.02; 0.02; 0.02; 0.02]; % Radius of bars [# bars x 1]
classK_test.Sradius = [0.01; 0.01]; % Radius of strings [# strings x 1]
classK_test.Nradius = 0.035*ones(size(classK_test.N,2),1); % Radius of node spheres [# nodes x 1]

%% Plotting configurations
% number_of_configurations = 3; % Number of configurations to plot
% tenseg_plot_configurations(Hist, classK_test, number_of_configurations)

%% Create animation
tenseg_animation(History,classK_test,[],[],[],[],10)

save Example_class_k_Dbar_Constraint.mat


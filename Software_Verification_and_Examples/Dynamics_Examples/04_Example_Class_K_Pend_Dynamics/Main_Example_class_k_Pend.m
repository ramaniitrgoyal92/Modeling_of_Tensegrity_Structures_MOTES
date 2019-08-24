% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 

clear all; clc; close all; warning off;

% Specify initial node positions for class k structure

%   - Note that you don't have to worry about converting class k nodes into
%     constrained class 1 physical and virtual nodes yourself.

n1 = [0 0 0]';
n2 = [1 0 0]';
N_simple = [n1 n2];

% Specify bar and string connectivity
C_b_in = [1 2];
C_s_in = [1 2];

% Convert above index notation into actual connectivity matrices
C_b = tenseg_ind2C(C_b_in,N_simple);
C_s = tenseg_ind2C(C_s_in,N_simple);

%%
% pinned_nodes
pinned_nodes = [1];

% Plot
tenseg_plot(N_simple,C_b,C_s);
%%

% Convert specified class k structure into a class 1 structure with constraints
[N_new,C_b_new,C_s_new,P,D,node_constraints] = tenseg_class_k_convert(N_simple,C_b,C_s,pinned_nodes);
% C_b_new = tenseg_ind2C(C_b_in_new,N_new);
% C_s_new = tenseg_ind2C(C_s_in_new,N_new);


% Print the final number of nodes
disp(['Converted Class K # of nodes:' num2str(size(N_new,2))])

% Print the generated node constraints
disp(['Class K Node constraints:'])
for i=1:numel(node_constraints)
    if numel(node_constraints{i})>1
        disp(['Coincident nodes: ' num2str(node_constraints{i})])
    end
end
%%
%Specify resting string lengths

% Here, we're setting every string rest length to 70% of its given length
S_0_percent = [(1:size(C_s_new,1))',1*ones(size(C_s_new,1),1)]; % percent of initial lengths

% This function converts those specified percentages into rest lengths
s_0 = tenseg_percent2s0(N_new,C_s_new,S_0_percent);


%%
% Add Velocity
V=zeros(3,length(N_new(1,:)));
% V(3,1:length(N(1,:))/3)=-2*ones(1,length(N(1,:))/3);
V(2,2) = -10;

% Add external force
W=zeros(3,length(N_new(1,:)));
W(2,2) = -10;



%% Simulation

% Create data structure of system BEFORE segmentation
classK_test.N = N_new;
classK_test.C_b = C_b_new;
classK_test.C_s = C_s_new;
classK_test.P = P;
classK_test.D = D;
classK_test.s_0 = s_0;
classK_test.tf = 5;
classK_test.dt = .001;
classK_test.video = 0;
classK_test.Nd0= V;
classK_test.W= W;

%%
% Perform simulation
[History,debug] = tenseg_sim_classkopen(classK_test);

%%
% Plot Node History
tenseg_plot_node(History,[1 2],[1 2 3])

% save Example_class_k_Pend_Dynamics.mat

%%
% Create animation
tenseg_animation(History,classK_test,[],[],[],[],10)




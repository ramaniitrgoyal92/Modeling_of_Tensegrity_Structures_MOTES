% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 

% EXAMPLE:
clear all; clc; close all;

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

% [N,C_b,C_s] = tenseg_prismplate(1,1,1);

% [N,C_b,C_s,info] = tenseg_skin(N,C_b,C_s,1);

prism.Bradius = 0.02*ones(size(C_b,1),1); % Radius of bars [# bars x 1]
prism.Sradius = 0.01*ones(size(C_s,1),1); % Radius of strings [# strings x 1]
prism.Nradius = 0.005*ones(size(N,2),1); % Radius of node spheres [# nodes x 1]

tenseg_plot(N,C_b,C_s,[],[],[],[],prism)


%%
% Divide strings into a number of segments
segments = 3;
[N,C_b,C_s,parents] = tenseg_string_segment(N,C_b,C_s,segments);

% Plot structure to see what string segmentation looks like
plot_handle = tenseg_plot(N,C_b,C_s);


%%
%Specify resting string lengths

% Here, we're setting every string rest length to 70% of its given length
S_0_percent = [(1:size(C_s,1))', 0.95*ones(size(C_s,1),1)]; % percent of initial lengths

% This function converts those specified percentages into rest lengths
s_0 = tenseg_percent2s0(N,C_s,S_0_percent,parents);

% W = 10*[zeros(2,size(N,2));[1 1 1 -1 -1 -1]];

%%
% Create data structure defining structure simulation task
prism.N = N;
% prism.Nd0 = ones(size(N));
prism.C_b = C_b;
prism.C_s = C_s;
prism.s_0 = s_0;
prism.c = 10;
% prism.W = W;

% These are optional inputs. If unspecified, default values will be loaded.
%    check tenseg_defaults.m to see/change these default values
prism.video = 0; % This turns off automatic animation output
prism.tf = 5;    % This sets the final simulation time to 3 sec 
prism.dt = 0.001; % This sets the simulation time step to 0.02 sec
prism.mb = 1;
prism.ms = 0.05;

% Perform Simulation
% History is a data structure containing simulation results
% sim_debug is a data structure containing a bunch of internally used variables
[History,sim_debug] = tenseg_sim_class1open(prism);


%% Create plots with 3D objects 
prism.Bradius = 0.02*ones(size(prism.C_b,1),1); % Radius of bars [# bars x 1]
prism.Sradius = 0.01*ones(size(prism.C_s,1),1); % Radius of strings [# strings x 1]
prism.Nradius = 0.02*ones(size(prism.N,2),1); % Radius of node spheres [# nodes x 1]

%% Create animation
tenseg_animation(History,prism,[],[],[],[],20)


%%
% Plot coordinate time histories for nodes 
tenseg_plot_node(History,[2 3 4],[1 3])    % plot x and y coordinates for nodes 2,3,4
                                           % not super useful because all
                                           % same colors right now...
save Example_class1_Prism_Dynamics.mat


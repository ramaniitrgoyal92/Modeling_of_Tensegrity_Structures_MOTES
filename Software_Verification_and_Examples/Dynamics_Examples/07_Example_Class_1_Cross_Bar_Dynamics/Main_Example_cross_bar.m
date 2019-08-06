
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 

% EXAMPLE:
clear all; clc; close all;

% Specify node positions
N = [-1 0 0; 0 -1 0; 1 0 0; 0 1 0]';

% Specify bar connectivity
Cb_in = [1 3; 2 4];
C_b = tenseg_ind2C(Cb_in,N);

% Specify string connectivity
Cs_in = [1 2; 2 3; 3 4;4 1];  % String one is node 1 to 2
C_s = tenseg_ind2C(Cs_in,N);

% Plot structure before string segmentation
tenseg_plot(N,C_b,C_s);

%%
% Divide strings into a number of segments
segments = 1;
[N,C_b,C_s,parents] = tenseg_string_segment(N,C_b,C_s,segments);

% Plot structure to see what string segmentation looks like
plot_handle = tenseg_plot(N,C_b,C_s);


%%
%Specify resting string lengths

% Here, we're setting every string rest length to 70% of its given length
S_0_percent = [(1:size(C_s,1))', [0.7*ones(size(C_s,1)/4,1);0.5*ones(size(C_s,1)/4,1);0.7*ones(size(C_s,1)/4,1);0.5*ones(size(C_s,1)/4,1)]]; % percent of initial lengths

% This function converts those specified percentages into rest lengths
s_0 = tenseg_percent2s0(N,C_s,S_0_percent,parents);


%%
% Create data structure defining structure simulation task
Cross_bar.N = N;
Cross_bar.C_b = C_b;
Cross_bar.C_s = C_s;
Cross_bar.s_0 = s_0;

% These are optional inputs. If unspecified, default values will be loaded.
%    check tenseg_defaults.m to see/change these default values
Cross_bar.video = 0; % This turns off automatic animation output
Cross_bar.tf = 5;    % This sets the final simulation time to 3 sec 
Cross_bar.dt = 0.001; % This sets the simulation time step to 0.02 sec
Cross_bar.c = 0;

% Perform Simulation
[History,sim_debug] = tenseg_sim_class1open(Cross_bar);

%%
% Create animation
tenseg_animation(History,Cross_bar,[],[],[],[],10)

figure()
tenseg_plot_node(History,[2 3 4],[1 2]) 

% save Example_class_k_X_Dynamics.mat

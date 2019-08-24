%% Example: String_Only_Dynamics

% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 

% EXAMPLE:
%    - Manual node matrix specification
%    - Manual connectivity matrix generation
%    - Barless structure
%    - Manual string rest length specification
%    - Time-varying external forces (see Wtest.m)
%    - Simulation

% Initialization
clear all; clc; close all; warning off;

%% Specify node positions

% Manually input node coordinates [x y z]'
n1 = [0 0 0]';
n2 = [1 0 0]';
n3 = [1 1 0]';
n4 = [0 1 0]';
% Build node matrix
N = [n1 n2 n3 n4]; 

%% Specify and build string connectivity matrix

% The first row [1 2] means string1 is from node 1 to 2, 
% The second row [2 3] means string2 is from node 2 to 3.
C_s_in = [1 2; 2 3; 3 4; 4 1; 2 4; 1 3]; 
% Connectivity matrix for strings
C_s = tenseg_ind2C(C_s_in,N); 

%% Plot structure

% Since there is no bars in this structure, C_b = []
tenseg_plot(N,[],C_s);

%% Specify resting string lengths

% This sets string 4 and 6 rest lengths to 90% and 100% of their given lengths
s_0 = tenseg_percent2s0(N,C_s,[4 0.9; 6 1]);

%% Create tensegrity simulation data structure

% Structure of net named as Snet
Snet.N = N;
Snet.C_s = C_s;
Snet.C_b = [];
Snet.ms = 1; % set point mass values
Snet.video = 0; % This turns off automatic animation output
% Snet.W = 'W_example'; % specify script that defines external forces at each time step
Snet.s_0 = s_0;
Snet.dt = 0.001; % ODE solver time step 
Snet.tf = 5; % Simulation time

% Perform simulation
[History,sim_debug] = tenseg_sim_class1open(Snet);

%% Plot History
% Plots simulation time histories of axes [1 3] (x and z) coordinate values
% for nodes [2 3 4]
tenseg_plot_node(History,[2 3 4],[1 3])

%% Create animation
tenseg_animation(History,Snet,[],[],[],[],10)



% save Example_String_Only_Dynamics.mat

%% ^_^ Welcome to Tensegrity Engineering Analysis Master(TEAM) software! ^_^ %%

%% SETUP file to be run only the first time

% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 

%% Add all necessary functions to MATLAB path
clear all; clc; close all
addpath 'Main' 'Test_Examples' 'Tensegrity_Examples' 'User_Guide' 'Videos' 'JOSS_Paper';
% Add subfolders of 'Main'
cd Main;
addpath 'Tensegrity_Statics_Main' 'Tensegrity_Dynamics_Main';
cd ..;
%% Open the User_Guide
cd User_Guide;
open('User_Guide_Tensegrity_Engineering_Analysis_Master_(TEAM).pdf');
cd ..;
disp('Welcome! Please follow the step-by-step instructions from the User Guide.');
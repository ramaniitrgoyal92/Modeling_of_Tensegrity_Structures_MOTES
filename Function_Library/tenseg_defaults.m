function default = tenseg_defaults(n,beta,alpha,sigma)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 
%
% default = TENSEG_DEFAULTS(n,beta,alpha,sigma) creates a structure
% containing all default variable specifications in the event that the
% variable is not explicitly specified by the user
%
% Inputs:
%	n: number of nodes
%	beta: number of bar members
%	alpha: number of string members
%	sigma: number of string point mass nodes
%
% Outputs:
%	default: structure containing default variable values

% Structure definitions
default.m = ones(1,beta); % bar masses
default.ms = 0.1*ones(1,sigma); % string node point masses
default.mgyro = zeros(1,beta); % gyro wheel masses
default.gyro_r = zeros(1,beta); % gyro wheel radii
default.gyro_h = zeros(1,beta); % gyro wheel heights
default.pinned_nodes = []; % indices of pinned nodes
default.P = []; % constraint matrix
default.D = []; % constraint matrix

% Initial conditions
default.Nd0 = zeros(3,n); %  node velocities

% Forces
default.W = zeros(3,n); % external node forces
default.k = 100*ones(1,alpha); % string stiffness coefficients
default.c = 0*ones(1,alpha); % string damping coefficients
default.gyro_omega = zeros(1,beta); % gyro angular rates

% Integration/Simulation
default.dt = 0.01; % integration time step
default.tf = 5; % integration final time
default.video = 1; % save video file toggle

end

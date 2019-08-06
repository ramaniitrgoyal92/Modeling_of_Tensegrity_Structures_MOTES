function [History,info] = tenseg_sim_class1open(sim_task)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 

% [History,info] = TENSEG_SIM_CLASS1OPEN(tenseg_in) performs simulation
% of a given CLASS 1, OPEN LOOP tensegrity simulation task. Specification
% of the simulation task is done with the 'sim_task' input data structure,
% described below. Outputs include a data structure containing node
% position time histories, and a data structure containing relevant
% internal values used in performing the simulation.
%
% Inputs:
%	sim_task: data structure describing simulation task
%		[].N: initial node positions
%		[].C_b: bar connectivity (empty if no bars)
%		[].C_s: string connectivity
%		[].s_0: string rest lengths
%
% Outputs:
%	History: data structure containing simulation results
%		[].Nhist: node positions for each time step
%		[].Ndhist: node velocities for each time step
%		[].gamma: string member force densities for each time step
%		[].t: time steps
%
% Example:
%	task.N = N;
%	task.C_b = C_b;
%	task.C_s = C_s;
%	task.s_0 = tenseg_percent2S0(N,C_s,0.8);
%	[hist,info] = tenseg_sim_class1open(task);


% Generate full sim task data structure from inputs
info = tenseg_struct_gen(sim_task);

% Get and log time span info
dt = info.sim.dt;
tf = info.sim.tf;
t = 0:dt:tf;
info.sim.tspan = t;

n = info.constants.n;


% ODE4 METHOD -------------------
% Reshape initial N, Nd matrices into vector format for ode function
XN0 = reshape(info.ICs.N0,3*n,1);
XNd0 = reshape(info.ICs.Nd0,3*n,1);
X0 = [XN0; XNd0];

% Store structure information to INPUTS to send to ode function
INPUTS.ICs = info.ICs;
INPUTS.tenseg_struct = info.sys_def;
INPUTS.forces = info.forces;
INPUTS.constants = info.constants;
INPUTS.sim = info.sim;

% Perform dynamics integration
Xhist = ode4(@tenseg_dyn_c1open_fnc,t,X0,INPUTS);
% Xhist = ode4C1(@tenseg_dyn_c1open_fnc,t,X0,INPUTS);

% Extract internally logged variables
logs = extractSignals(@tenseg_dyn_c1open_fnc,{'bar_len','gamma'});

% Reshape and store simulation results
Nhist = reshape(Xhist(:,1:3*n)',3,n,numel(t));
Ndhist = reshape(Xhist(:,3*n+1:end)',3,n,numel(t));
History.Nhist = Nhist;
History.Ndhist = Ndhist;
History.bar_len = logs.bar_len;
History.gamma = logs.gamma;
History.t = t;


% Make animation video
if info.sim.video
	tenseg_animation(History,info.sys_def)
end

end
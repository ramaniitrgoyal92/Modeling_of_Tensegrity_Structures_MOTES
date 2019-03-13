function [History,info] = tenseg_sim_classkopen(sim_task)

% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 

% [History,info] = TENSEG_SIM_CLASSKOPEN(tenseg_in) performs simulation
% of a given CLASS K, OPEN LOOP tensegrity simulation task. Specification
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
%		[].P: constraint matrix
%       [].D: constraint matrix
%
% Outputs:
%	History: data structure containing simulation results
%		[].Nhist: node positions for each time step
%		[].Ndhist: node velocities for each time step
%		[].gamma: string member force densities for each time step
%		[].t: time steps
%
% Example:
%	[N,C_b,C_s,P] = tenseg_class_k_convert(N,C_b,C_s);
%	task.N = N;
%	task.C_b = C_b;
%	task.C_s = C_s;
%	task.s_0 = tenseg_percent2S0(N,C_s,0.8);
%	task.P = P;
%	[hist,info] = tenseg_sim_class1open(task);


% Generate full sim task data structure from inputs
info = tenseg_struct_gen(sim_task);


% Get and log time span info
dt = info.sim.dt;
tf = info.sim.tf;
t = 0:dt:tf;
info.sim.tspan = t;

n = info.constants.n;
P = info.sys_def.P;
D = info.sys_def.D;

%     For Reduced Order
    [U,SIGMA,V] = svd(P);
    rank_SIGMA = rank(SIGMA);
    U2 = U(:,rank_SIGMA+1:end);
    V1 = V(:,1:rank_SIGMA);
    Sigma1 = SIGMA(1:rank_SIGMA,1:rank_SIGMA); 	
    ETA1 = D*V1*Sigma1^-1;
    
    
% ODE4 METHOD -------------------
% Resphape initial condition matrices into vectors to pass into ODE4
% XN0 = reshape(info.ICs.N0,3*n,1);
% XNd0 = reshape(info.ICs.Nd0,3*n,1);
% X0 = [XN0; XNd0];

ETA20 = info.ICs.N0*U2;
ETA2d0 = info.ICs.Nd0*U2;
XETA20 = reshape(ETA20,[],1);
XETA2d0 = reshape(ETA2d0,[],1);
X0 = [XETA20;XETA2d0];

% Create INPUTS structure to pass variables into ODE4 integration
INPUTS.ICs = info.ICs;
INPUTS.tenseg_struct = info.sys_def;
INPUTS.forces = info.forces;
INPUTS.constants = info.constants;
INPUTS.sim = info.sim;

% Perform simulation
% Xhist = ode4(@tenseg_dyn_ckopen_fnc,t,X0,INPUTS);
Xhist = ode4Ck(@tenseg_dyn_ckopen_fnc,t,X0,INPUTS);

% Extract internally logged variables and label them
logs = extractSignals(@tenseg_dyn_ckopen_fnc,{'bar_len','gamma'});

ETA2hist = reshape(Xhist(:,1:size(Xhist,2)/2)',3,[],numel(t));
ETA2dhist = reshape(Xhist(:,size(Xhist,2)/2+1:end)',3,[],numel(t));

% Reshape simulation results back into matrix form
for i = 1:numel(t)
    Nhist(:,:,i) = [ETA1 ETA2hist(:,:,i)]*U';
    Ndhist(:,:,i) = [zeros(size(ETA1)) ETA2dhist(:,:,i)]*U';
end

% Log simulation results
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
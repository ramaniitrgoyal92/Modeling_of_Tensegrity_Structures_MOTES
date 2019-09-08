function [full_struct] = tenseg_struct_gen(input_struct)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 
%
% [tenseg_struct] = TENSEG_STRUCT_GEN(tenseg_struct) fully populates all
% fields of a tensegrity simulation task data structure given an initial
% data structure containing, at a minimum, the required fields. Default
% values are loaded and added as fields where necessary.
%
% Inputs:
%	input_struct: tensegrity simulation task specification data structure
%		REQUIRED FIELDS:
%			[].N: initial node positions
%			[].C_b: bar connectivity
%			[].C_s: string connectivity
%			[].s_0: resting string lengths
%		OPTIONAL FIELDS (default values specified in 'tenseg_defaults'):
%			[].P: constraint matrix
%			[].D: constraint matrix
%			[].Nd0: initial node velocities
%			[].m: bar masses
%			[].ms: string node point masses
%			[].mgyro: gyro wheel masses
%			[].gyro_r: gyro wheel radii
%			[].gyro_h: gyro wheel heights
%			[].tf: simulation time duration
%			[].dt: simulation time step
%			[].W: external node forces
%			[].k: string stiffness coefficients
%			[].gyro_omega: gyro angular rates
%
% Outputs:
%	full_struct: summary of internal variables used
%		[].input: input variables given (see input field descriptions)
%		[].sys_def: tensegrity system definition variables
%			[].[].C_sb: string connections to bars
%			[].[].C_ss: string connections to string nodes
%			[].[].C_bb: non-zero columns of C_b (old C_b matrix)
%			[].[].C_nb: satisfies R_b = N*C_nb^T*C_r'
%			[].[].C_s: overall string connectivity
%			[].[].C_b: overall bar connectivity
%			[].[].m: bar masses
%			[].[].ms: string node point masses
%			[].[].mgyro: gyro wheel masses
%			[].[].gyro_r: gyro wheel radii
%			[].[].gyro_h: gyro wheel heights
%			[].[].C_r: Bar center of mass 'connectivity'
%			[].[].len_hat: initial bar lengths
%			[].[].Jt_hat: transaxial bar+gyro moments of inertia (/len?)
%			[].[].Ja_hat: axial bar+gyro moments of inertia
%		[].constants: various constants throughout the simulation
%			[].[].n: number of nodes
%			[].[].beta: number of bars
%			[].[].alpha: number of string members
%			[].[].sigma: number of string point mass nodes
%			[].[].M: mass matrix used in dynamics
%			[].[].Minv: inverse of M
%		[].forces: force-related variables
%			[].[].W: external node forces
%			[].[].k: string stiffness coefficients
%			[].[].s_0: resting string lengths
%			[].[].gyro_omega: gyro angular rates
%		[].ICs: initial conditions
%			[].[].N0: initial node positions
%			[].[].Nd0: initial node velocities
%		[].sim: integration specification variables
%			[].[].dt: simulation time step
%			[].[].tf: simulation time duration


% Save original input fields
full_struct.input = input_struct;


% Read in mandatory fields ------------------------------------------------
N = input_struct.N;
C_b = input_struct.C_b;
C_s = input_struct.C_s;


% Get basic structure values
n = size(N,2); % number of nodes
beta = size(C_b,1); % number of bars
alpha = size(C_s,1); % number of strings

% Check if we have string point mass nodes
sigma = n-2*beta;
% if sigma==0,
% 	% Reorder nodes to canonical form
% 	[N,C_b,C_s] = tenseg_N_reorder(N,C_b,C_s);
% end

% Load all default input values
default = tenseg_defaults(n,beta,alpha,sigma);



% Initialize initial conditions struct ------------------------------------
% Mandatory user inputs
ICs.N0 = N; % Initial node positions

% Default to specified default values if unspecified
ICs = tenseg_field_init('Nd0',input_struct,ICs,default,0); % Initial node velocities
Nd0 = ICs.Nd0;
full_struct.ICs = ICs;

% Check for Correct Initial Velocity 

if ~isempty(C_b)
check_velocity = diag((N*C_b')'*(Nd0*C_b'));

if norm(check_velocity) >= 1e-5
    disp('Incorrect Initial Velocities');
    return;
end
end

% Initialize structure definitions struct ----------------------------------
% Mandatory user inputs
sys_def.N = N;
sys_def.C_s = C_s; % String connectivity
sys_def.C_b = C_b; % Bar Connectivity

sys_def = tenseg_convert_Cmats(sys_def);
C_r = sys_def.C_r; % Bar center of mass 'connectivity'
C_nb = sys_def.C_nb;
C_ns = sys_def.C_ns;
C_bb = sys_def.C_bb;

% Default to specified default values if unspecified
sys_def = tenseg_field_init('P',input_struct,sys_def,default,0); % constraint matrix
sys_def = tenseg_field_init('D',input_struct,sys_def,default,0); % constraint matrix
sys_def = tenseg_field_init('m',input_struct,sys_def,default,0); % bar masses
sys_def = tenseg_field_init('ms',input_struct,sys_def,default,0); % string node masses
sys_def = tenseg_field_init('mgyro',input_struct,sys_def,default,0); % gyro masses
sys_def = tenseg_field_init('gyro_r',input_struct,sys_def,default,0); % gyro wheel radii
sys_def = tenseg_field_init('gyro_h',input_struct,sys_def,default,0); % gyro wheel heights

% If only one point mass value is specified, assume all point masses have
% that value
if numel(sys_def.ms)==1
	sys_def.ms = sys_def.ms*ones(1,sigma);
end

% If only one bar mass value is specified, assume all bar masses have
% that value
if numel(sys_def.m)==1
	sys_def.m = sys_def.m*ones(1,beta);
end


% Check if structure has bar members
if beta
	sys_def.len_hat = sqrt(diag(diag((N*C_b')'*(N*C_b')))); % Bar lengths diagonal matrix
	len = diag(sys_def.len_hat)';
else
    sys_def.len_hat = 0;
	len = 0;
end
	
% Moments of inertia ---------------
gyro_r = sys_def.gyro_r;
gyro_h = sys_def.gyro_h;
k_val = 1./(12*len.^2).*(3*gyro_r.^2 + gyro_h.^2);

mb = sys_def.m;
mgyro = sys_def.mgyro;
Jt = mb/12 + k_val.*mgyro;
sys_def.Jt_hat = diag(Jt);

Ja = 1/2*mgyro./len.*gyro_r.^2;
sys_def.Ja_hat = diag(Ja);

full_struct.sys_def = sys_def;



% Initialize constants struct ---------------------------------------------
constants.n = n;
constants.beta = beta;
constants.alpha = alpha;
constants.sigma = sigma;

m_tot_hat = diag(sys_def.m + sys_def.mgyro);
m_s_hat = diag(sys_def.ms);
Jt_hat = sys_def.Jt_hat;
Ja_hat = sys_def.Ja_hat;

M = [C_nb'*(C_bb'*Jt_hat*C_bb + C_r'*m_tot_hat*C_r), C_ns'*m_s_hat];
Minv = M^-1;

constants.M = M;
constants.Minv = Minv;

full_struct.constants = constants;



% Initialize simulation constants struct ----------------------------------
sim = [];
sim = tenseg_field_init('dt',input_struct,sim,default,0); % Integration time step value
sim = tenseg_field_init('tf',input_struct,sim,default,0); % Integration final time
sim = tenseg_field_init('video',input_struct,sim,default,0); % Create video afterwards

full_struct.sim = sim;



% Initialize force variables struct ---------------------------------------
forces = [];
forces = tenseg_field_init('W',input_struct,forces,default,0); % External forces acting on each node
forces = tenseg_field_init('k',input_struct,forces,default,0); % String stiffness coefficients
forces = tenseg_field_init('c',input_struct,forces,default,0); % String damping coefficients
forces.s_0 = input_struct.s_0; % Resting string lengths
forces = tenseg_field_init('gyro_omega',input_struct,forces,default,0); % gyro wheel angular rates

if numel(forces.k)==1
	forces.k = forces.k*ones(1,alpha);
end
if numel(forces.c)==1
	forces.c = forces.c*ones(1,alpha);
end

full_struct.forces = forces;
end




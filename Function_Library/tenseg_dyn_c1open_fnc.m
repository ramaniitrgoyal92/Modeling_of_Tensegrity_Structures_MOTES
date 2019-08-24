function [ Xd, logging] = tenseg_dyn_c1open_fnc( t,X,INPUTS )
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 

% [ Xd, logging] = TENSEG_DYN_C1OPEN_FNC( t,X,INPUTS ) is used in computing
% the dynamical behavior of an open-loop class 1 tensegrity structure.
% Essentially, X_dot = [Nd; Ndd] is computed from X = [N, Nd], and t is
% included for use with ode4, a fixed time-step RK4 integrator.
% Structure-specific definitions are provided through INPUTS, and
% functionality is included that allows for the logging of internally
% computed variables throughout a simulation, such as bar length and string
% force densities. At a high level, this internal logging works by storing
% the internal variables in a persistent internal variable. If the function
% is then called without any inputs, the logged internal values are
% returned.
%
% Inputs:
%	t: current time (for use with ode4)
%	X: current N and N_dot values
%		stacked column-wise into a vector: X = [N; Nd];
%	INPUTS: structure containing tensegrity structure-specific information
%		below fields required from TENSEG_SIM_CLASS1OPEN.m
%		[].constants
%		[].sim
%		[].forces
%		[].tenseg_struct
%
% Outputs:
%	Xd: N_dot and N_dotdot values based on current N and N_dot. Again,
%		stacked column-wise into a vector Xd = [Nd; Ndd]
%	logging: used for outputting internally logged variables


% Used to log internal variables. Remain in workspace through multiple
% function calls
persistent i logging_temp


% Prep for logging internal variables -------------------------------------
% Labels for internally-logged data
bar_len_log = 1;
gamma_log = 2;

% Check if this is the start of a simulation and initialize logging arrays
if nargout == 1 && (isempty(i) || i > 2) && t==0
	i = 1;
	
	% Read in some inputs to initialize size of internal logging arrays
	constants = INPUTS.constants;
	sim = INPUTS.sim;
	
	% Initialize persistent logging arrays
	logging_temp{bar_len_log} = zeros(constants.beta,numel(sim.tspan));
	logging_temp{gamma_log} = zeros(constants.alpha,numel(sim.tspan));
end


% If simulation (one argument output), read inputs and compute Xd
if nargout == 1
	
	%disp(t)
	
	idx = 0.25 * (i + 3); % used for checking if we're on a reported time-step
	ind_use = floor((i-2)/4)+1; % index to use when reading in time-varying info
	if ind_use==0
		ind_use=1;
	end
	dt = INPUTS.sim.dt;
	
	% Read inputs
	forces = INPUTS.forces;
	tenseg_struct = INPUTS.tenseg_struct;
	constants = INPUTS.constants;

	Wgiven = forces.W;
	C_b = tenseg_struct.C_b;
	C_s = tenseg_struct.C_s;
    len_hat = tenseg_struct.len_hat;    
	m = tenseg_struct.m;
	Minv = constants.Minv;
	n = constants.n;
	
	% Reshape from vector to matrix form
	XN = X(1:numel(X)/2);
	XNd = X(numel(X)/2+1:end);
	N = reshape(XN,3,n);
	Nd = reshape(XNd,3,n);

    %         Bar Length Correction - Use only if necessary - Needs to be verified
    if ~isempty(C_b)
%         [N,Nd] = bar_length_correction_Class1(N,Nd,C_b,len_hat);
    end
    
	% Get current external forces
	if ischar(Wgiven)
		run(Wgiven) % if W is a string, run that script
	elseif size(Wgiven,3)==1
		W = Wgiven; % W can be constant
	else
		W = Wgiven(:,:,ind_use); % or W can be time-varying
    end
    
    

	% Compute bar lengths
	% -- needs to be done at each time step so numerical errors are at least
	%    self-consistent
	if ~isempty(C_b) % for string-only structure, C_b should be empty
		B = N*C_b';
	else
		B = 0;
    end
    Current_len_hat = sqrt(diag(diag(B'*B)));


	% Computing gamma -----------------------------------------------------
	S = N*C_s';
	Sd = Nd*C_s';

	% Initialize diagonal gamma matrix
	gamma_hat = zeros(size(C_s,1));

	% Get current resting string lengths
	if ischar(forces.s_0)
		run(forces.s_0); % s_0 can be a script name
	elseif size(forces.s_0,2)==1
		S_0_hat = diag(forces.s_0); % s_0 can be constant
	else
		S_0_hat = diag(forces.s_0(:,ind_use)); % or s_0 can be time-varying
	end

	% Go through each string and figure out its force density
	for j=1:size(C_s,1)
		if norm(S(:,j))>S_0_hat(j,j) % makes sure tension in string is 0 if it's shorter than rest length
			gamma_hat(j,j) = forces.k(1,j)*(1-S_0_hat(j,j)/norm(S(:,j)));
			damping_term = forces.c(1,j)*(S(:,j)'*Sd(:,j))/norm(S(:,j))^2;
			gamma_hat(j,j) = gamma_hat(j,j) + damping_term;            
                    
        if gamma_hat(j,j) < 0
            gamma_hat(j,j) = 0;
        end 
        
        end    
	end
	gamma_hat = sparse(gamma_hat);


	
	% LOG values of interest -- NOTE: specific to ode4 --------------------
	if ~rem(idx, 1)
        % disp(t)
		logging_temp{gamma_log}(:,idx) = diag(gamma_hat);
		logging_temp{bar_len_log}(:,idx) = diag(Current_len_hat);
	end
	
	

	% CLASS 1 Dynamics ----------------------------------------------------
	if ~isempty(C_b) % If string-only strcture, set Bd to 0
		Bd = Nd*C_b';
	else
		Bd = 0;
	end

	Ja_hat = tenseg_struct.Ja_hat;
	gyro_omega_hat = diag(forces.gyro_omega);
	C_bb = tenseg_struct.C_bb;

	% Make Bcross matrix
	beta = constants.beta;
	Bcross = zeros(3,3*beta);
	for ii=1:beta
		Bcross(:,(ii-1)*3+1:(ii-1)*3+3) = skew(B(:,ii));
	end

	% Make Bd_hat matrix
	Bd_hat = zeros(3*beta,beta);
	for ii=1:beta
		Bd_hat((ii-1)*3+1:(ii-1)*3+3,ii) = Bd(:,ii);
	end

	Q = Bcross*Bd_hat;

	sigma = constants.sigma;
	Wsgyro_current = W + [Q*Ja_hat*gyro_omega_hat*Current_len_hat^-2*C_bb, zeros(3,sigma)];
	F = (Wsgyro_current-N*C_s'*gamma_hat*C_s);

	Jt_hat = tenseg_struct.Jt_hat;
	
	% If string-only structure, set lambda_hat to 0
	if ~isempty(C_b)
		lambda_hat = -Jt_hat*Current_len_hat^-2*diag(diag(Bd'*Bd))-1/2*Current_len_hat^-2*diag(diag(B'*F*C_b'));
	else
		lambda_hat = 0;
	end
		
	C_sb = tenseg_struct.C_sb;
	C_ss = tenseg_struct.C_ss;
	C_nb = tenseg_struct.C_nb;
	K = [C_s'*gamma_hat*C_sb - C_nb'*C_bb'*lambda_hat*C_bb, C_s'*gamma_hat*C_ss];

	% Calculate node accelerations
	Ndd = (Wsgyro_current-N*K)*Minv;

	% Reshape Nd, Ndd for output
	XNdd = reshape(Ndd,3*n,1);
	Xd = [XNd; XNdd];

	i=i+1;
end

% If this function is called with more than 1 output argument, it returns
% the internally logged data
if nargout > 1
	Xd = [];
	logging = logging_temp;
end

end


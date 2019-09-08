function [ Xd, logging] = tenseg_dyn_ckopen_fnc( t,X,INPUTS )
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 
%
% [ Xd, logging] = TENSEG_DYN_CKOPEN_FNC( t,X,INPUTS ) is used in computing
% the dynamical behavior of an open-loop class k tensegrity structure.
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


% Temporary switch between analytical and numerical Omega solver.
% Analytical seems to only be working for Class 2...?
numerical_omega_solver = 0;

% Prep for logging internal variables -------------------------------------
persistent i logging_temp

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

% If simulation, read inputs and compute Xd
if nargout == 1
	% disp(t)
	
	idx = 0.25 * (i + 3); % used for checking if we're on a reported time-step
	ind_use = floor(i/4)+1; % index to use when reading in time-varying info
	
	% Read inputs
	forces = INPUTS.forces;
	tenseg_struct = INPUTS.tenseg_struct;
	constants = INPUTS.constants;

	Wgiven = forces.W;
	C_b = tenseg_struct.C_b;
	C_s = tenseg_struct.C_s;
	P = tenseg_struct.P;
	D = tenseg_struct.D;
	len_hat = tenseg_struct.len_hat;    
	constraints = size(P,2);
	m = tenseg_struct.m;
	Minv = constants.Minv;
	M = constants.M;
	n = constants.n;
	
    
%     For Reduced Order
    [U,SIGMA,V] = svd(P);
    rank_SIGMA = rank(SIGMA);
    U1 = U(:,1:rank_SIGMA);
    U2 = U(:,rank_SIGMA+1:end);
    V1 = V(:,1:rank_SIGMA);
    Sigma1 = SIGMA(1:rank_SIGMA,1:rank_SIGMA); 	
    ETA1 = D*V1*Sigma1^-1;
   
    
	% Reshape from vector to matrix form
    
	XETA2 = X(1:numel(X)/2);
	XETA2_dot = X(numel(X)/2+1:end);
	ETA2 = reshape(XETA2,3,[]);
	ETA2_dot = reshape(XETA2_dot,3,[]);
    
    N = [ETA1 ETA2]*U';
    Nd = [zeros(size(ETA1)) ETA2_dot]*U';
    
    
    %     Bar Length Correction - Use only if necessary - Needs to be verified
    [N,Nd] = bar_length_correction_ClassK(N,Nd,C_b,P,D,len_hat);

    ETA2 = N*U2;
    
	% Compute bar lengths
	% -- needs to be done at each time step so numerical errors are at least
	%    self-consistent
	if ~isempty(C_b) % Check if there are bars in the system
		B = N*C_b';
	else
		B = 0;
	end
	Current_len_hat = sqrt(diag(diag(B'*B)));


	% Get current external forces
	if ischar(Wgiven)
		run(Wgiven) % if W is a string, run that script
	elseif size(Wgiven,3)==1
		W = Wgiven; % W can be constant
	else
		W = Wgiven(:,:,ind_use); % or W can be time-varying
	end


	% Compute gamma -----------------------------------------------------
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
		logging_temp{gamma_log}(:,idx) = diag(gamma_hat);
		logging_temp{bar_len_log}(:,idx) = diag(Current_len_hat);
	end
	
	
	% CLASS K Dynamics ----------------------------------------------------
	
	% Compute B_dot if the structure has bars
	if ~isempty(C_b)
		Bd = Nd*C_b';
	else
		Bd = 0;
	end

	Ja_hat = tenseg_struct.Ja_hat;
	gyro_omega_hat = diag(forces.gyro_omega);
	C_bb = tenseg_struct.C_bb;

	% Make Bcross matrix (for gyro dynamics)
	beta = constants.beta;
	Bcross = zeros(3,3*beta);
	for ii=1:beta
		Bcross(:,(ii-1)*3+1:(ii-1)*3+3) = skew(B(:,ii));
	end

	% Make Bd_hat, Q matrices (for gyro dynamics)
	Bd_hat = zeros(3*beta,beta);
	for ii=1:beta
		Bd_hat((ii-1)*3+1:(ii-1)*3+3,ii) = Bd(:,ii);
	end
	Q = Bcross*Bd_hat;
 
    
	% NUMERICAL METHOD -- REMOVE ONCE CONFIDENT IN ANALYTICAL SOLUTION
	% Simultaneously solve for Ndd and Omega (constraint lagrange multipliers)
	if numerical_omega_solver
		cvx_begin quiet
			variable Omega(3,constraints)

			sigma = constants.sigma;
			Wsgyro_current = W + [Q*Ja_hat*gyro_omega_hat*Current_len_hat^-2*C_bb, zeros(3,sigma)] + Omega*P';
			F = (Wsgyro_current-N*C_s'*gamma_hat*C_s);

			Jt_hat = tenseg_struct.Jt_hat;
			if ~isempty(C_b)
				lambda_hat = -Jt_hat*Current_len_hat^-2*diag(diag(Bd'*Bd))-1/2*Current_len_hat^-2*diag(diag(B'*F*C_b'));
			else
				lambda_hat = 0;
			end

			C_sb = tenseg_struct.C_sb;
			C_ss = tenseg_struct.C_ss;
			C_nb = tenseg_struct.C_nb;
			K = [C_s'*gamma_hat*C_sb - C_nb'*C_bb'*lambda_hat*C_bb, C_s'*gamma_hat*C_ss];

			Ndd = (Wsgyro_current-N*K)*Minv;

			minimize(norm(Omega))
			subject to
				abs(Ndd*P) <= 10*eps*zeros(size(Ndd*P))
		cvx_end
    else
        
		% ANALYTICAL METHOD1 - Full Order
		
		m_hat = diag(m);

		bC=P'*C_b';
		bD=C_b*Minv*P;
		bE=P'*Minv*P;
		E=eye(3);
		bA=-S*gamma_hat*C_s*Minv*P+B*...
			diag(diag(.5*Current_len_hat^(-2)*B'*(S*gamma_hat*C_s-W)*C_b'-1/12*Current_len_hat^(-2)*m_hat*(Bd'*Bd)))...
			*C_b*Minv*P+W*Minv*P;
		Lag_Mat=0;
		for j=1:size(C_b,1)
			Lag_Mat=Lag_Mat+1/(2*Current_len_hat(j,j)^2)*kron(bC(:,j)',kron(B(:,j),(B(:,j)*bD(j,:))'));
		end
		Lag_Mat=Lag_Mat-[kron(bE,E(1,:));kron(bE,E(2,:));kron(bE,E(3,:))];

		OMEGA_LAGRANGE_Vec=Lag_Mat\[bA(1,:)';bA(2,:)';bA(3,:)'];
		OMEGA_LAGRANGE=reshape(OMEGA_LAGRANGE_Vec,3,numel(OMEGA_LAGRANGE_Vec)/3);
        Omega = OMEGA_LAGRANGE;
        
	% ANALYTICAL METHOD2 - Reduced Order
    
% 		m_hat = diag(m);
% 
% 		bC=P'*C_b';
% 		bD=C_b*Minv*U1;
% 		bE=P'*Minv*U1;
% 		E=eye(3);
% 		bA=-S*gamma_hat*C_s*Minv*U1+B*...
% 			diag(diag(.5*Current_len_hat^(-2)*B'*(S*gamma_hat*C_s-W)*C_b'-1/12*Current_len_hat^(-2)*m_hat*(Bd'*Bd)))...
% 			*C_b*Minv*U1+W*Minv*U1;
% 		Lag_Mat=0;
% 		for j=1:size(C_b,1)
% 			Lag_Mat=Lag_Mat+1/(2*Current_len_hat(j,j)^2)*kron(bC(:,j)',kron(B(:,j),(B(:,j)*bD(j,:))'));
% 		end
% 		Lag_Mat=Lag_Mat-[kron(bE,E(1,:));kron(bE,E(2,:));kron(bE,E(3,:))];
% 
% 		OMEGA_LAGRANGE_Vec=Lag_Mat\[bA(1,:)';bA(2,:)';bA(3,:)'];
% 		OMEGA_LAGRANGE=reshape(OMEGA_LAGRANGE_Vec,3,numel(OMEGA_LAGRANGE_Vec)/3);
%         Omega = OMEGA_LAGRANGE;
    


		% Using lagrange multipliers, solve for Ndd from dynamical equations
		C_sb = tenseg_struct.C_sb;
		C_ss = tenseg_struct.C_ss;
		C_nb = tenseg_struct.C_nb;
		sigma = constants.sigma;
		Jt_hat = tenseg_struct.Jt_hat;

           
        
		Wsgyro_current = W + [Q*Ja_hat*gyro_omega_hat*Current_len_hat^-2*C_bb, zeros(3,sigma)] + Omega*P';
		F = (Wsgyro_current-N*C_s'*gamma_hat*C_s);

		if ~isempty(C_b)
			lambda_hat = -Jt_hat*Current_len_hat^-2*diag(diag(Bd'*Bd))-1/2*Current_len_hat^-2*diag(diag(B'*F*C_b'));
		else
			lambda_hat = 0;
		end

		K = [C_s'*gamma_hat*C_sb - C_nb'*C_bb'*lambda_hat*C_bb, C_s'*gamma_hat*C_ss];
        M_bar = U2'*M*U2;
        K_bar = U2'*K*U2;
        W_bar = (W + [Q*Ja_hat*gyro_omega_hat*Current_len_hat^-2*C_bb, zeros(3,sigma)])*U2-ETA1*U1'*K*U2;
        
        ETA2_dotdot = (W_bar-ETA2*K_bar)/M_bar;
%         Ndd = [zeros(size(ETA1)) ETA2_dotdot]*U';

% 		Ndd = (Wsgyro_current-N*K)*Minv;
	end
	
	% Reshape Nd, Ndd for output
	XETA2_dotdot = reshape(ETA2_dotdot,[],1);
	Xd = [XETA2_dot; XETA2_dotdot];

	i=i+1;
end

% If this function is called with more than 1 output argument, it returns
% the internally logged data
if nargout > 1
	Xd = [];
	logging = logging_temp;
end


function tenseg_struct = tenseg_convert_Cmats(tenseg_struct)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 

% tenseg_struct = TENSEG_CONVERT_CMATS(tenseg_struct) converts given
% 'simple' connectivity matrices (C_b, C_s) into the the connectivity
% matrices used in the full string-to-string dynamics (C_bb, C_ns, C_nb,
% C_r)
%
% Inputs:
%	tenseg_struct: structure containing tensegrity system definitions
%		[].N
%		[].C_b
%		[].C_s
%
% Outputs:
%	tenseg_struct: modified structure containing the newly computed
%		connectivity matrices

N = tenseg_struct.N;
C_b = tenseg_struct.C_b;
C_s = tenseg_struct.C_s;


% Get basic structure values
n = size(N,2); % number of nodes
beta = size(C_b,1); % number of bars


% Check if we have string point mass nodes
sigma = n-2*beta;


tenseg_struct.N_b = N(:,1:2*beta);
tenseg_struct.N_s = N(:,2*beta+1:end);


% Convert to separate matrices for dynamics (also calc of consts M,Minv)
C_sb = C_s(:,1:2*beta); % string connections to bars
tenseg_struct.C_sb = C_sb;
C_ss = C_s(:,2*beta+1:end); % string connections to string nodes
tenseg_struct.C_ss = C_ss;

C_bb = C_b(:,1:2*beta); % non-zero columns of C_b (old C_b matrix)
tenseg_struct.C_bb = C_bb;

C_r = 1/2*abs(C_b(:,1:2*beta)); % Bar center of mass 'connectivity'
tenseg_struct.C_r = C_r;

C_nb = [eye(2*beta), zeros(2*beta,sigma)]; % satisfies R_b = N*C_nb^T*C_r'
tenseg_struct.C_nb = C_nb;
C_ns = [zeros(sigma,2*beta), eye(sigma)]; % satisfies R_s = N*C_ns'
tenseg_struct.C_ns = C_ns;
	
end
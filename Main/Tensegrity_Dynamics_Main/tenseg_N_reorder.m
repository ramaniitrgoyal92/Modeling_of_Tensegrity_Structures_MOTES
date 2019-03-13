function [N_out,C_b_out,C_s_out] = tenseg_N_reorder(N,C_b,C_s)

% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 

% COMPATIBILITY WITH OTHER FUNCTIONS NEEDS TO BE TESTED 
% [N_out,C_b_out,C_s_out] = TENSEG_N_REORDER(N,C_b,C_s) Rearranges nodes
% such that C_b simplifies to [eye(beta), -eye(beta)] for a structure with
% beta bars.
%
% Inputs:
%	N: initial node matrix (3 x n array for n nodes)
%	C_b: initial bar connectivity matrix (beta x n array for beta bars)
%	C_s: initial string connectivity matrix (alpha x n array for alpha)
%		strings)
%
% Outputs:
%	N_out: rearranged node matrix
%	Cb_out: rearranged bar connectivity matrix
%	Cs_out: rearranged string connectivity matrix

beta = size(C_b,1); % Number of bars

% Initialize sorting vector for bar nodes
new_N_order = zeros(1,2*beta);

for i=1:beta % Go through all rows in Cb
	
	% Find indices for current bar
	side1_ind = find(C_b(i,:)==1);
	side2_ind = find(C_b(i,:)==-1);
	
	% Store indices in order required to get desired C_b form
	new_N_order(i) = side1_ind;
	new_N_order(i+beta) = side2_ind;
end

% Rearrange N, C_b, and C_s
N_out = N(:,new_N_order);
C_b_out = C_b(:,new_N_order);
C_s_out = C_s(:,new_N_order);

end
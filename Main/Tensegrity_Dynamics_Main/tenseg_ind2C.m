function [ C_mat ] = tenseg_ind2C( C_ind,Nodes )
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 


% [ C_mat ] = TENSEG_IND2C( C_ind, Nodes ) creates a connectivity matrix
% from input index notation array and node matrix.
%
% Inputs:
%	C_ind: index connectivity array (m x 2 array for m members)
%	Nodes: node matrix (3 x n array for n nodes)
%
% Outputs:
%	C_mat: connectivity matrix (m x n matrix satisfying M = N*C')
%
% Example: If given four nodes (N is a 3x4 matrix), and you want one bar to
%	be the vector from node 1 to 3, and another bar from node 2 to 4, input
%	C_ind would be: C_ind = [1 3;
%						     2 4]
%
% C_b = tenseg_ind2C([1 3; 2 4], N);


nmembers = size(C_ind,1); % Number of members being created
n = size(Nodes,2); % Number of nodes in the structure

% Initialize connectivity matrix
C_mat = zeros(n,nmembers);


for i=1:nmembers % Go through each member
	
	% Get indices for start and end points
	side1=C_ind(i,1);
	side2=C_ind(i,2);
	
	% Put -1 at index for start point, 1 at index for end point
	C_mat(side1, i) = -1;
	C_mat(side2, i) = 1;
end
C_mat = C_mat';

end

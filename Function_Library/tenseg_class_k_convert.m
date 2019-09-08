function [N_new,C_b_new,C_s_new,P,D,node_constraints] = tenseg_class_k_convert(N,C_b,C_s,pinned_nodes)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 
%
% [N_new,C_b_new,C_s_new,P,node_constraints] =
% TENSEG_CLASS_K_CONVERT(N,C_b,C_s,pinned_nodes) converts a given class k
% tensegrity structure into a corresponding class 1 structure with
% constrained physical and virtual nodes
% 
% Inputs:
%	N: Initial (class k) node matrix
%	C_b: Initial bar connectivity matrix
%	C_s: Initial string connectivity matrix
%	pinned_nodes (optional): indices of nodes that are pinned in place
%
% Outputs:
%	N_new: Converted node matrix with physical and virtual nodes
%	C_b_new: Converted bar connectivity matrix ''
%	C_s_new: Converted string connectivity matrix ''
%	P: Constraint matrix satisfying N*P = D
%	node_constraints: Cell array in which each cell gives constrained node
%		indices corresponding to the initial nodes

if nargin < 4
	pinned_nodes = []; % default pinned nodes
end

% Get class of each node, and total number of virtual nodes
node_class = sum(abs(C_b),1);
number_of_virtual_nodes = sum(node_class)-sum(node_class>0);

n = numel(node_class);

% Initialize matrix of virtual nodes
N_virtual = zeros(3,number_of_virtual_nodes);

% Go through all nodes. If a node has one or more coincident virtual nodes,
% copy that node position into N_virtual. Also create node_constraints
% along the way, which keeps track of which final nodes are constrained to
% be conincident
virtual_node_ind = 1;
for i=1:numel(node_class) % go through each given node
	node_constraints{i}(1) = i;
	for j=1:node_class(i)-1 % go through any virtual nodes for node i
		N_virtual(:,virtual_node_ind) = N(:,i);
		node_constraints{i}(1+j) = n + virtual_node_ind; % log new index
		virtual_node_ind = virtual_node_ind+1;
	end
end
N_new = [N N_virtual]; % combine physical and virtual node matrices

% Go through each row of given bar connectivity matrix and move entries for
% class k nodes to the appropriate virtual node indices
C_b_new = C_b;
for i=1:numel(node_constraints)
	for j=1:numel(node_constraints{i})-1
		node_ind = node_constraints{i}(j);
		if node_ind <= size(C_b,2)
			bar_ind = find(abs(C_b(:,node_ind))~=0);
			for k=1:numel(bar_ind)-1
				C_b_new(bar_ind(k+1),node_constraints{i}(k+1)) = C_b_new(bar_ind(k+1),node_ind);
				C_b_new(bar_ind(k+1),node_ind) = 0;
			end
		end
	end
end

% Pad C_s with zeros to match dimension of new N matrix
C_s_new = [C_s, zeros(size(C_s,1),number_of_virtual_nodes)];


% Generate P constraint matrix
P = [];
for i=1:numel(node_constraints)
	if numel(node_constraints{i})>1
		constraints_to_add_more = nchoosek(node_constraints{i},2);
        constraints_to_add = constraints_to_add_more(1:size(node_constraints{i},2)-1,:);
%         constraints_to_add = constraints_to_add_more;
		P_add = zeros(size(N_new,2),size(constraints_to_add,1));
		for j=1:size(constraints_to_add,1)
			P_add(constraints_to_add(j,1),j) = 1;
			P_add(constraints_to_add(j,2),j) = -1;
		end
		P = [P P_add];
	end
end


% Rearrange stuff to put into form: N = [Nb Ns]
string_nodes_ind = [node_class==0, logical(zeros(1,number_of_virtual_nodes))];
N_s = N_new(:,string_nodes_ind);
N_b = N_new(:,~string_nodes_ind);
N_new = [N_b N_s];

C_b_s = C_b_new(:,string_nodes_ind);
C_b_b = C_b_new(:,~string_nodes_ind);
C_b_new = [C_b_b C_b_s];

C_s_s = C_s_new(:,string_nodes_ind);
C_s_b = C_s_new(:,~string_nodes_ind);
C_s_new = [C_s_b C_s_s];



% Generate D matrix
D = zeros(3,size(P,2));
for i=1:numel(pinned_nodes)
	D_add = N(:,pinned_nodes(i));
	D = [D D_add];
end

% Generate constraint matrix for pinned nodes
P_pinned = [];
for i=1:numel(pinned_nodes)
	P_add = zeros(size(N_new,2),1);
	P_add(pinned_nodes(i)) = 1;
	P_pinned = [P_pinned P_add];
end

P = [P P_pinned];

if ~isempty(P)
	P_s = P(string_nodes_ind,:);
	P_b = P(~string_nodes_ind,:);
	P = [P_b; P_s];
end

end

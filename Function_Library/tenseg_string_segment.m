function [N_new,C_b_new,C_s_new,parents,ms,seg_k,seg_c] = tenseg_string_segment(N,C_b,C_s,segments,string_mass_assign,stiffness,damping)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 
%
% [N_new,C_b_new,C_s_new,parents,ms] = TENSEG_STRING_SEGMENT(N,C_b,C_s,segments,string_mass_assign)
% splits each string member of a given tensegrity structure into a
% specified number of segments.
%
% Inputs:
%	N: initial node matrix (3 x n for n nodes)
%	C_b: initial bar connectivity matrix (beta x n for beta bars)
%	C_s: initial string connectivity matrix (alpha x n for alpha strings)
%	segments: number of segments that each string member is divided into
%		if segments = [k], all strings will be divided into k segments
%		if segments = [ind k], the ind-th string will be divided
%			into k segments, and all else will be left alone
%	string_mass_assign (optional): Either a m x 2 matrix specifying string
%		masses (first column giving string index, second giving
%		mass), or a single value specifying mass for all strings
%	stiffness (optional): Either a m x 2 matrix specifying string
%		stiffness values (first column giving string index, second giving
%		stiffness), or a single value specifying stiffness for all strings
%	damping (optional): Either a m x 2 matrix specifying string
%		damping values (first column giving string index, second giving
%		damping), or a single value specifying damping for all strings
%
% Outputs:
%	N_new: updated node matrix containing string point mass nodes
%	C_b_new: updated bar connectivity matrix
%	C_s_new: updated string connectivity matrix
%	parents: vector containing original string indices associated with
%		each generated string segment
%	ms: point mass values determined from given string masses


if nargin<5
	string_mass_assign = [];
end
if nargin<6
	stiffness = [];
end
if nargin<7
	damping = [];
end


alpha = size(C_s,1);

% If there's a single number in 'segments', every string will be split into
% that many segments
if numel(segments)==1
	assign_segments = segments*ones(1,alpha);
else
	assign_segments = ones(1,alpha);
	assign_segments(segments(:,1)) = segments(:,2);
end

% If segments is 1, we can save some time and just leave the input
% unchanged.
if segments==1
	N_new = N;
	C_b_new = C_b;
	C_s_new = C_s;
	parents = 1:size(C_s,1);
else
	S = N*C_s'; % Get all string vectors

	N_s = []; % initialize string point mass matrix
	C_s_s = []; % initialize string-to-string connectivity matrix
	parents = []; % initialize parent info logging

	segment_ind = 1;
	
	% Loop through each string to segment
	for j=1:alpha
		
		% Get starting node for new string segment
		node_segment_start_ind = find(C_s(j,:)==-1);
		node_segment_start = N(:,node_segment_start_ind);

		C_s_s(segment_ind, node_segment_start_ind) = -1;

		parents = [parents, j];
		
		if assign_segments(j)>1
			
			% Get vector of current segment
			segment_vector = S(:,j)/assign_segments(j);
			
			% Loop for each segment on current string
			for i=1:assign_segments(j)-1
				% Log parent string index for current segment
				parents = [parents, j];

				% Get new string point mass node
				new_node = node_segment_start + segment_vector;
				N_s = [N_s, new_node];
				C_s_s(segment_ind, size(N,2)+size(N_s,2)) = 1;

				% Set new node as start node for next segment
				node_segment_start = new_node;

				segment_ind = segment_ind+1;
				C_s_s(segment_ind, size(N,2)+size(N_s,2)) = -1;

				if i==assign_segments(j)-1
					original_string_end_ind = find(C_s(j,:)==1);
					C_s_s(segment_ind, original_string_end_ind) = 1;
					segment_ind = segment_ind+1;
				end
			end
		else
			node_segment_start_ind = find(C_s(j,:)==-1);
			C_s_s(segment_ind, node_segment_start_ind) = -1;
			
			original_string_end_ind = find(C_s(j,:)==1);
			C_s_s(segment_ind, original_string_end_ind) = 1;
			segment_ind = segment_ind+1;
		end
	end

	N_new = [N, N_s]; % full node matrix
	C_s_new = C_s_s; % full string connectivity matrix
	C_b_new = [C_b, zeros(size(C_b,1),(size(C_s_new,2)-size(C_b,2)))]; % bar connectivity using N
end


% Compute point mass values based on given string masses
% If single mass value given, convert that into more general form
if numel(string_mass_assign)==1
	string_mass_assign = [(1:size(C_s,1))' string_mass_assign*ones(size(C_s,1),1)]; 
end

% Divide string mass evenly between string point masses
ms = [];
for i=1:size(string_mass_assign,1)
	current_ind = string_mass_assign(i,1);
	number_of_nodes = assign_segments(i)-1;
	point_mass_val = string_mass_assign(i,2)/number_of_nodes;
	ms_add = point_mass_val*ones(1,number_of_nodes);
	ms = [ms ms_add];
end


% Compute segment stiffness values based on given parent string stiffnesses
% If single stiffness value given, convert that into more general form
if numel(stiffness)==1
	stiffness = [(1:size(C_s,1))' stiffness*ones(size(C_s,1),1)]; 
end

% Calculate stiffness terms for each segment such that parent stiffness is maintained
seg_k = [];
for i=1:size(stiffness,1)
	current_ind = stiffness(i,1);
	number_of_segments = assign_segments(i);
	seg_k_val = stiffness(i,2)*number_of_segments;
	seg_k_add = seg_k_val*ones(1,number_of_segments);
	seg_k = [seg_k seg_k_add];
end


% Compute segment damping values based on given parent string damping vals
% If single damping value given, convert that into more general form
if numel(damping)==1
	damping = [(1:size(C_s,1))' damping*ones(size(C_s,1),1)]; 
end

% Calculate damping terms for each segment such that parent damping is maintained
seg_c = [];
for i=1:size(damping,1)
	current_ind = damping(i,1);
	number_of_segments = assign_segments(i);
	seg_c_val = damping(i,2)/number_of_segments;
	seg_c_add = seg_c_val*ones(1,number_of_segments);
	seg_c = [seg_c seg_c_add];
end


end

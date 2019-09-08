function [s_0] = tenseg_percent2s0(N,C_s,percents,parents)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 
%
% s_0 = TENSEG_PERCENT2S0(N,C_s,percents,parents) creates a vector
% containing all resting string lengths for a set of string segments based
% on specified rest length percentages specified for the "parent" strings.
%
% Inputs:
%	N: node matrix
%	C_s: string connectivity matrix (after segmentation)
%	percents: m x 2 matrix for m parent strings that have resting lengths
%		that differ from their specified initial lenghts. The first column
%		gives the index of the mth parent string, and the second column
%		gives the percent of the initial length that defines the rest
%		length. If values are being specified for every string, indices can
%		be left out. If all strings have same percent, can be single value.
%	parents (optional): REQUIRED if working with segmented strings. Vector
%		containing original string indices associated with each string
%		segment (from 'tenseg_string_segment'). Not needed otherwise.
%
% Outputs:
%	s_0: vector containing rest lengths of all string segments
%
% Example:
%	s_0 = tenseg_percent2s0(N,C_s,[1 0.8; 2 0.8; 3 0.7]);


% Each string is its own parent if not otherwise specified
if nargin==3 || isempty(parents)
	parents = 1:size(C_s,1);
end

% If single percent value given, assign it to every string member
if numel(percents)==1
	percents = percents*ones(numel(unique(parents)),1);
end

% If only percents given for all parent strings, generate indices
if size(percents,2)==1 && numel(percents)==numel(unique(parents))
	percents = [(1:numel(unique(parents)))' percents];
end


% Get all initial string lengths
s_0 = sqrt(diag((N*C_s')'*(N*C_s')));


% Compute ith string rest length
for i=1:size(percents,1)
	
	% Get parent string index
	p_ind = percents(i,1);
	
	% Get percentage by which to scale rest length
	scale = percents(i,2);
	
	% Scale all string members associated with parent string index
	s_0(parents==p_ind) = scale*s_0(parents==p_ind);
end

end

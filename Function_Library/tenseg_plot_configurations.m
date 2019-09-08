function tenseg_plot_configurations(History,tenseg_struct,...
    n_frames,highlight_nodes,view_vec)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 
%
% tenseg_plot_configurations(History,tenseg_struct,n_frames,highlight_nodes,view_vec)
% Creates plots for uniformly time-spaced configurations, the number of
% configurations is given by input n_frames
%
% Inputs:
%	History: simulation results from tenseg_sim_[] functions. Mainly just
%		using History.Nhist, which is a 3xnxt array containing node
%		positions for all n nodes in the structure at all t time steps in
%		the simulation   
%	tenseg_struct: tensegrity structure definition data structure
%		containing bar and string connectivity matrices
%	n_frames: Number of frames to plot
%	highlight_nodes (optional): node(s) to highlight when plotting (vector
%		containing node numbers)
%	view_vec (optional): plotting view direction (see view())

if nargin < 4
	highlight_nodes = [];
end
if nargin < 5
	view_vec = [];
end

Nhist = History.Nhist;
C_b = tenseg_struct.C_b;
C_s = tenseg_struct.C_s;

%% Collecting data for 3D plots
R3Ddata = [];
if isfield(tenseg_struct,'Bradius')
    R3Ddata.Bradius = tenseg_struct.Bradius;
end
if isfield(tenseg_struct,'Sradius')
    R3Ddata.Sradius = tenseg_struct.Sradius;
end
if isfield(tenseg_struct,'Nradius')
    R3Ddata.Nradius = tenseg_struct.Nradius;
end

%%
% Set axes for animation
[axis_vec, view_vec_derived] = tenseg_axisview(Nhist, R3Ddata);

% if no view vector defined, use the one derived from Nhist
if isempty(view_vec)
	view_vec = view_vec_derived;
end

% selecting time frames
frames = floor(linspace(1,size(History.t,2),n_frames));

% go through all time frames, plot, adjust axes
for i=1:n_frames
	N = Nhist(:,:,frames(i));
    ititle = ['t = ', num2str(History.t(frames(i)))];
	tenseg_plot(N,C_b,C_s,[],highlight_nodes,view_vec,ititle,R3Ddata);
	axis equal
	axis(axis_vec)

	drawnow
end

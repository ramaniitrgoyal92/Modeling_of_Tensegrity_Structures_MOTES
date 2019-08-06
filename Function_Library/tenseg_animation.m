function [] = tenseg_animation(History,tenseg_struct,filename,highlight_nodes,view_vec,realtime,frame_skip)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 

% [] =
% tenseg_animation(History,tenseg_struct,filename,highlight_nodes,view_vec,realtime)
% creates an animation video file for a given tensegrity sim history and
% the topology of the structure
%
% Inputs:
%	History: simulation results from tenseg_sim_[] functions. Mainly just
%		using History.Nhist, which is a 3xnxt array containing node
%		positions for all n nodes in the structure at all t time steps in
%		the simulation   
%	tenseg_struct: tensegrity structure definition data structure
%		containing bar and string connectivity matrices
%	filename (optional): path/filename for video file
%	highlight_nodes (optional): node(s) to highlight when plotting (vector
%		containing node numbers)
%	view_vec (optional): plotting view direction (see view())
%	realtime (optional): force animation video file to be in real-time
%		(1sec = 1sec). Default: 0 (No)
%   frame_skip (optional): number of frames to skip
%
% Example:
%	[hist,info] = tenseg_sim_class1open(prism);
%	tenseg_animation(hist,prism)


Nhist = History.Nhist; % Node matrix history
C_b = tenseg_struct.C_b; % Connectivity matrices
C_s = tenseg_struct.C_s;

% Handle optional inputs
if nargin < 3 || isempty(filename)
	filename = 'test'; % default file name
end
if nargin < 4
	highlight_nodes = []; % default nodes to highlight
end
if nargin < 5
	view_vec = []; % default view vector
end
if nargin < 6
	realtime = 0; % default realtime video framerate option
end
if nargin < 7
    frame_skip = 1;
end


% Initialize video save file
vid = VideoWriter(filename);

% Figure out framerate and number of steps to skip to make realtime
if realtime
	% Determine true framerate based on t_span
	true_framerate = numel(History.t)/History.t(end);
	
	% Limit framerate to 30 fps
	if true_framerate >= 30
		vid.FrameRate = 30;
		
		% Determine # of frames to skip based on framerate
		frame_skip = round(true_framerate/vid.FrameRate);
	else
		vid.FrameRate = round(true_framerate);
		frame_skip = 1;
	end
else
	% If not realtime, just print every frame
end

vid.Quality = 100;
open(vid);


% Create figure in which to plot each frame
%fig = figure();
fig = figure('units','normalized','position',[0 0 0.8 0.8]);

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
% Determine view/axes based on Node matrix history
[axis_vec, view_vec_derived] = tenseg_axisview(Nhist,R3Ddata);

% If no view vector defined, use the one derived from Nhist
if isempty(view_vec)
	view_vec = view_vec_derived;
end
	

% Go through all time steps, plot, adjust axes, and write to video file
for i=1:frame_skip:size(Nhist,3)
	N = Nhist(:,:,i);
	fig = tenseg_plot(N,C_b,C_s,fig,highlight_nodes,view_vec,[],R3Ddata);
	axis equal
	axis(axis_vec)
	
	% Get frame width/height and write frame (including axes) to file
	%position = get(gcf,'Position'); % do this to include axes in video
	%writeVideo(vid,getframe(gcf,[0 0 position(3:4)]));

	drawnow
	SelectEntireFig = getframe(fig);
	writeVideo(vid,SelectEntireFig);    
	clf
end
close(vid);

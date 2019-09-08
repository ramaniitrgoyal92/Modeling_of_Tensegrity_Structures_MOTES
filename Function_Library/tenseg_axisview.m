function [axis_vec, view_vec] = tenseg_axisview(Nhist, R3Ddata)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 
%
% [axis_vec, view_vec] = TENSEG_AXISVIEW(Nhist) finds axes bounds and an
% appropriate view vector for visualizing data from a provided node matrix
% or node history array.
%
% Inputs:
%	Nhist: node matrix or node history array (3 x n x t)
%   R3Ddata (optional): Structure with the radius of objects for 3D plots
%        [].Bradius: Radius of bars [# bars x 1]
%        [].Sradius: Radius of strings [# strings x 1]
%        [].Nradius: Radius of node spheres [# nodes x 1]
%
% Outputs:
%	axis_vec: min and max axes values
%	view_vec: viewing direction for visualization

if nargin < 2
    R3Ddata = [];
end

% Get min and max X,Y,Z coordinates contained in Nhist 
min_x = min(min(Nhist(1,:,:)));
max_x = max(max(Nhist(1,:,:)));
min_y = min(min(Nhist(2,:,:)));
max_y = max(max(Nhist(2,:,:)));
min_z = min(min(Nhist(3,:,:)));
max_z = max(max(Nhist(3,:,:)));

%%
maxradius = 0;
if ~isempty(R3Ddata) % Extending axis to account for radii of objects
    if isfield(R3Ddata,'Bradius')
        maxradius = max(maxradius,max(R3Ddata.Bradius));
    end
    if isfield(R3Ddata,'Sradius')
        maxradius = max(maxradius,max(R3Ddata.Sradius));
    end
    if isfield(R3Ddata,'Nradius')
        maxradius = max(maxradius,max(R3Ddata.Nradius));
    end   
    min_x = min_x - maxradius; max_x = max_x + maxradius;
    min_y = min_y - maxradius; max_y = max_y + maxradius;
    min_z = min_z - maxradius; max_z = max_z + maxradius;
end

%%
% Get difference between min and max values for each axis
diff_x = max_x-min_x;
diff_y = max_y-min_y;
diff_z = max_z-min_z;

% If values vary in three dimensions, set orthogonal view with
% corresponding min/max axis values
if diff_x>eps+2*maxradius && diff_y>eps+2*maxradius && diff_z>eps+2*maxradius
	view_vec_derived = [1 1 1];
	axis_vec = [min_x max_x min_y max_y min_z max_z];
	
% If values only vary in two dimensions, set view for 2D plotting
else
	if diff_z < eps+2*maxradius
		view_vec_derived = [0 0 1];
		axis_vec = [min_x max_x min_y max_y];
	elseif diff_y < eps+2*maxradius
		view_vec_derived = [0 1 0];
		axis_vec = [min_x max_x 0 1 min_z max_z];
	elseif diff_x < eps+2*maxradius
		view_vec_derived = [1 0 0];
		axis_vec = [0 1 min_y max_y min_z max_z];
	end
end
view_vec = view_vec_derived;

end

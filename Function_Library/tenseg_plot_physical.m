function [fig_out] = tenseg_plot_physical( N,C_b,C_s,fig_handle,highlight_nodes,view_vec, PlotTitle, R3Ddata,Yld_n_Bkl,BarWidth_in,StringWidth_in)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 
%
% This function added thickness display apart from tenseg_plot.m
% Thickness of bars and strings are proportional to the force densities in
% them.
%
% [fig_out] = TENSEG_PLOT( N,C_b,C_s,fig_handle,highlight_nodes,view_vec )
% creates a rough visualization figure for a given tensegrity structure
%
% Inputs:
%	N: node matrix (3 x n array for n nodes)
%	C_b (optional): bar connectivity matrix (beta x n array for beta bars)
%	C_s: string connectivity matrix (alpha x n array for alpha strings)
%
%   LABEL: Gives the information on bar failure mode; 11/23/18
%
%	fig_handle (optional): figure in which to plot
%	highlight_nodes (optional): node(s) to highlight when plotting (vector
%		containing node numbers)
%	view_vec (optional): plotting view direction (see view())
%   PlotTitle (optional): Title of the plot, default: ''
%   R3Ddata (optional): Structure with the radius of objects for 3D plots
%        [].Bradius: Radius of bars [# bars x 1]
%        [].Sradius: Radius of strings [# strings x 1]
%        [].Nradius: Radius of node spheres [# nodes x 1]
%   Yld_n_Bkl: Yielding or Buckling label
%   BarWidth_in: line width of bars
%   StringWidth_in: line width of strings
% Outputs:
%	fig_out: figure handle in which plot was made
%
% Example:
%	tenseg_plot_physical( N,C_b,C_s,fig_handle,highlight_nodes,view_vec, PlotTitle, R3Ddata,Yld_n_Bkl,BarWidth_in,StringWidth_in)

% Handle optional arguments

%% Object size options (for line plots)
% BarWidth = BarWidth_in; % Width of bar lines
% StringWidth = StringWidth_in ; % Width of string lines
NodeSize = 20; % Size of node marker

%% Labeling options
% Write labels? (1: show, 0: suppress)
LabelNodes = 1;
LabelStrings = 1;
LabelBars = 1;

FontBars = 15; % Font of bar labels
FontStrings = 10; % Font of string labels
FontNodes = 18; % Font of node labels

FractionDistance = 0.005; % Distance between object and label (relative to overall axis size)

%% 3D plot options
nsurfpatches = 6; % Number of surface patches in 3D plots
BarSurfColor = [0.2, 0.2, 0.6];
StringSurfColor = [0.9, 0.1, 0.1];
LightAmbientStrength = 0.7;% [0,1], 0.3 is Matlab's default

%%
switch nargin
    case 5
        view_vec = [];
    case 4
        view_vec = [];
        highlight_nodes = [];
    case 3
        view_vec = [];
        highlight_nodes = [];
        fig_handle = [];
end

if nargin < 7 % Empty Plot Title
    PlotTitle = '';
end
if nargin < 8 % No 3D object data
    R3Ddata = [];
end

% Get min and max X,Y,Z coordinates contained in Nhist
min_x = min(N(1,:,:)); max_x = max(N(1,:,:));
min_y = min(N(2,:,:)); max_y = max(N(2,:,:));
min_z = min(N(3,:,:)); max_z = max(N(3,:,:));

%%
if ~isempty(R3Ddata) % Extending axis to account for radii of objects
    maxradius = 0;
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

% Get difference between min and max values for each axis
diff_x = max_x-min_x;
diff_y = max_y-min_y;
diff_z = max_z-min_z;

% Get distance for labels
dist_x = FractionDistance * diff_x;
dist_y = FractionDistance * diff_y;
dist_z = FractionDistance * diff_z;

%% Open specified figure or create new one
if isempty(fig_handle)
    fig_out = figure;
else
    fig_out = figure(fig_handle);
end


%% Plot bar member vectors
if ~isempty(C_b)
    B = N*C_b';
    bar_start_nodes = zeros(3,size(B,2));
    bar_end_nodes = zeros(3,size(B,2));
    for j = 1:size(B,2)
        bar_start_nodes(:,j) = N(:,C_b(j,:)==-1);
        bar_end_nodes(:,j) = N(:,C_b(j,:)==1);
    end
    
    if ~isempty(R3Ddata) && isfield(R3Ddata,'Bradius') % 3D plot
        for j = 1:size(C_b,1)
            [LatFace, UpFace, DwFace] = PlotCylinderObject(bar_start_nodes(:,j),bar_end_nodes(:,j),...
                R3Ddata.Bradius(j),nsurfpatches); % Lateral, upper, and down surface of cylinder representing a bar
            
            surf(LatFace.x, LatFace.y, LatFace.z, 'MeshStyle','row','FaceColor',BarSurfColor, ...
                'FaceLighting','gouraud', 'AmbientStrength',LightAmbientStrength);
            hold on
            surf(UpFace.x, UpFace.y, UpFace.z, 'MeshStyle','row','FaceColor',BarSurfColor, ...
                'FaceLighting','gouraud', 'AmbientStrength',LightAmbientStrength);
            hold on
            surf(DwFace.x, DwFace.y, DwFace.z, 'MeshStyle','row','FaceColor',BarSurfColor, ...
                'FaceLighting','gouraud', 'AmbientStrength',LightAmbientStrength);
            hold on
        end
    else % Normal line plot
        
        for jj=1:size(B,2)
            quiver3(bar_start_nodes(1,jj),bar_start_nodes(2,jj),bar_start_nodes(3,jj),B(1,jj),B(2,jj),B(3,jj),'black.','Autoscale','off','LineWidth',BarWidth_in(jj))
            hold on
        end
        
    end
    % Write bar labels
    for i = 1:size(B,2)
        text(bar_start_nodes(1,i) + 0.5*B(1,i) + dist_x, ...
            bar_start_nodes(2,i) + 0.5*B(2,i)+0.1 + dist_y, ...
            bar_start_nodes(3,i) + 0.5*B(3,i) + dist_z,...
            num2str(i), 'FontSize',FontBars, 'Color', 'k')
        
    end
    
    Bkl=find(Yld_n_Bkl==1); Yld=find(Yld_n_Bkl==0);
    
    for j=1:length(Bkl)
        text(bar_start_nodes(1,Bkl(j)) + 0.5*B(1,Bkl(j))+0.1 + dist_x, ...
            bar_start_nodes(2,Bkl(j)) + 0.5*B(2,Bkl(j))+0.1 + dist_y, ...
            bar_start_nodes(3,Bkl(j)) + 0.5*B(3,Bkl(j)) + dist_z,...
            'Buckling', 'FontSize',FontBars, 'Color', 'k')
    end
    
    for k=1:length(Yld)
        text(bar_start_nodes(1,Yld(k)) + 0.5*B(1,Yld(k))+0.1 + dist_x, ...
            bar_start_nodes(2,Yld(k)) + 0.5*B(2,Yld(k))+0.1 + dist_y, ...
            bar_start_nodes(3,Yld(k)) + 0.5*B(3,Yld(k)) + dist_z,...
            'Yielding', 'FontSize',FontBars, 'Color', 'k')
    end
    
end

%% Plot string member vectors
if ~isempty(C_s)
    S = N*C_s';
    string_start_nodes = zeros(3,size(S,2));
    string_end_nodes = zeros(3,size(S,2));
    for j = 1:size(S,2)
        string_start_nodes(:,j) = N(:,C_s(j,:)==-1);
        string_end_nodes(:,j) = N(:,C_s(j,:)==1);
    end
    
    if ~isempty(R3Ddata) && isfield(R3Ddata,'Sradius') % 3D plot
        for j = 1:size(C_s,1)
            [LatFace, UpFace, DwFace] = PlotCylinderObject(string_start_nodes(:,j),string_end_nodes(:,j),...
                R3Ddata.Sradius(j),nsurfpatches); % Lateral, upper, and down surface of cylinder representing a string
            
            surf(LatFace.x, LatFace.y, LatFace.z, 'MeshStyle','row','FaceColor',StringSurfColor, ...
                'FaceLighting','gouraud', 'AmbientStrength',LightAmbientStrength);
            hold on
            surf(UpFace.x, UpFace.y, UpFace.z, 'MeshStyle','row','FaceColor',StringSurfColor, ...
                'FaceLighting','gouraud', 'AmbientStrength',LightAmbientStrength);
            hold on
            surf(DwFace.x, DwFace.y, DwFace.z, 'MeshStyle','row','FaceColor',StringSurfColor, ...
                'FaceLighting','gouraud', 'AmbientStrength',LightAmbientStrength);
            hold on
        end
    else % Normal line plot
        
        for kk=1:size(S,2)
            quiver3(string_start_nodes(1,kk),string_start_nodes(2,kk),string_start_nodes(3,kk),S(1,kk),S(2,kk),S(3,kk),'red.','Autoscale','off','LineWidth',StringWidth_in(kk));
            hold on
        end
    end
    
    % Write string labels
    if LabelStrings == 1
        for i = 1:size(S,2)
            text(string_start_nodes(1,i) + 0.45*S(1,i) + dist_x, ...
                string_start_nodes(2,i) + 0.45*S(2,i) + dist_y, ...
                string_start_nodes(3,i) + 0.45*S(3,i) + dist_z,...
                num2str(i), 'FontSize',FontStrings, 'Color', 'r')
        end
    end
end

%% Plot node markers
if ~isempty(R3Ddata) && isfield(R3Ddata,'Nradius') % 3D plot
    for i = 1:size(N,2)
        [xnod,ynod,znod] = sphere(nsurfpatches);
        xnod = R3Ddata.Nradius(i)*xnod + N(1,i);
        ynod = R3Ddata.Nradius(i)*ynod + N(2,i);
        znod = R3Ddata.Nradius(i)*znod + N(3,i);
        surf(xnod,ynod,znod,...
            'FaceColor','k','EdgeColor','none','FaceLighting','gouraud',...
            'AmbientStrength',LightAmbientStrength);
        hold on
    end
else % Normal plot
    plot3(N(1,:),N(2,:),N(3,:),'black.','MarkerSize',NodeSize)
end

% Write node labels
if LabelNodes == 1
    for i = 1:size(N,2)
        text(N(1,i) + dist_x, N(2,i) + dist_y, N(3,i) + dist_z,...
            num2str(i), 'FontSize', FontNodes, 'Color', 'b')
    end
end

% Highlight specified nodes if applicable
for j=1:numel(highlight_nodes)
    node_index = highlight_nodes(j);
    plot3(N(1,node_index),N(2,node_index),N(3,node_index),'rd','MarkerSize',8,'MarkerFaceColor','red')
end

% Modify plot display
grid off
axis equal
if isempty(view_vec)
    [~, view_vec_derived] = tenseg_axisview(N,R3Ddata);
    view_vec = view_vec_derived;
end
view(view_vec)

xlabel('x')
ylabel('y')
zlabel('z')
title(PlotTitle)
set(gca,'fontsize', 15,'linewidth',1.15)
set(gca,'ticklength',1.2*get(gca,'ticklength'))

if ~isempty(R3Ddata) % 3D object changes
    camlight('headlight');
    camproj('perspective');
end


end

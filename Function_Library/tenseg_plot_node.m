function [] = tenseg_plot_node(History,node_ind,axes)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 

% tenseg_plot_node(History,node_ind,axes,line_style) plots time histories
% of specified axes coordinate values for a given node from a given
% simulation time history. Line style can optionally be specified.
%
% Inputs:
%	History: structure containing [].Nhist obtained from tenseg_sim funciton
%	node_ind: index of node to plot
%	axes (optional): coordinate axes to plot (defaults to [1 2 3])
%
% Example:
%	tenseg_plot_node(hist,1,1:2)

figure()
% Handle optional argument inputs
switch nargin
	case 2
		axes = 1:3;
end

if isempty(axes)
	axes = 1:3;
end

number_of_plots = numel(axes);

% Create subplots for each axis
for i=1:number_of_plots
	subplot(number_of_plots,1,i)
	plot(History.t, squeeze(History.Nhist(axes(i),node_ind,:))','LineWidth',2)
	hold on
	if axes(i)==1
		ylabel('$x$','Interpreter','latex')
	elseif axes(i)==2
		ylabel('$y$','Interpreter','latex')
	elseif axes(i)==3
		ylabel('$z$','Interpreter','latex')
    end
    set(gca,'fontsize', 15,'linewidth',1.15)
    set(gca,'ticklength',1.2*get(gca,'ticklength'))
    for i = 1:size(node_ind,2)
    legendInfo{i} = ['Node' '$' num2str(node_ind(1,i)) '$'];
end
legend(legendInfo,'Interpreter','latex')
xlabel('Time','Interpreter','latex')
end



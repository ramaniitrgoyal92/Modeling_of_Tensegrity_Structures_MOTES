%% Double Pendulum Simulation Analytical Solution
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 

clc;clear all;close all;
% Initial Conditions
m1 = 1; m2 = 1;
L1 = 1; L2 = 1;
theta1 = 1*pi/4;
theta2 = 1*pi/4;
t = linspace(0,5,50001);

%% Solving ODE of a double pendulum
[T,Y] = ode45(@(t,x) double_pendulum_dynamics(t,x,m1,L1,m2,L2),t, [theta1,theta2,0,0]);

%% Calculating joint coordinates for animation purposes
x = [L1*sin(Y(:,1)),L1*sin(Y(:,1))+L2*sin(Y(:,2))];
y = [-L1*cos(Y(:,1)),-L1*cos(Y(:,1))-L2*cos(Y(:,2))];

% %% Plot the results
% % Convert radians to degrees
% L = L1;
% angle = Y(:,1:2)*180/pi;
% % Set up first frame
% figure('Color', 'white')
% % subplot(2,1,1)
% plot(T, angle,'LineWidth',2)
% hh1(1) = line(T(1),angle(1,1),'Marker','.','MarkerSize',20,'Color','b');
% hh1(2) = line(T(1),angle(1,2),'Marker','.','MarkerSize',20,'Color','r');
% xlabel('Time(sec)')
% ylabel('Angle(deg)')
% 
% figure
% % subplot(2,1,2)
% hh2 = plot([0,x(1,1);x(1,1),x(1,2)], [0,y(1,1);y(1,1),y(1,2)],'.-','MarkerSize',20,'LineWidth',2);
% axis equal
% axis([-2*L 2*L -2*L 2*L])
% ht = title(sprintf('Time: %0.2f sec',T(1)));
% % Get figure size
% pos = get(gcf,'Position');
% width = pos(3);
% height = pos(4);
% 
% % Preallocate data (for storing frame data)
% mov = zeros(height,width,1,length(T));%, 'uint8');
% % Loop through by changing XData and YData
% for id = 1:length(T)
%    % Update graphics data. This is more efficient than recreating plots.
%    set(hh1(1),'XData',T(id),'YData',angle(id,1))
%    set(hh1(2),'XData',T(id),'YData',angle(id,2))
%    set(hh2(1),'XData',[0,x(id,1)],'YData',[0,y(id,1)])
%    set(hh2(2),'XData',x(id,:),'YData',y(id,:))
%    set(ht,'String',sprintf('Time: %0.2f sec',T(id)))
% 
%    % output gif animation    
%      frame=getframe(gcf);  
%      im=frame2im(frame);% make gif file
%      [I,map]=rgb2ind(im,20);
%      if id ==1
%          imwrite(I,map,'double_pendulum_animation.gif','gif','Loopcount',inf,'DelayTime',0);% create when 1st time run
%      else
%        imwrite(I,map,'double_pendulum_animation.gif','gif','WriteMode','append','DelayTime',0);
%      end
% end

%%
figure('Color', 'white')
subplot(2,1,1)
plot(T, zeros(size(x(:,1))), T, x(:,1),T, x(:,2),'LineWidth',2)
xlabel('Time','Interpreter','latex')
ylabel('$x$','Interpreter','latex')
set(gca,'fontsize', 15,'linewidth',1.15)
set(gca,'ticklength',1.2*get(gca,'ticklength'))
legend('Node 1','Node 2','Node 3','Interpreter','latex')
subplot(2,1,2)
plot(T, zeros(size(y(:,1))), T, y(:,1),T, y(:,2),'LineWidth',2)
xlabel('Time','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
legend('Node 1','Node 2','Node 3','Interpreter','latex')
set(gca,'fontsize', 15,'linewidth',1.15)
set(gca,'ticklength',1.2*get(gca,'ticklength'))

%%
save Analytical_Solution.mat

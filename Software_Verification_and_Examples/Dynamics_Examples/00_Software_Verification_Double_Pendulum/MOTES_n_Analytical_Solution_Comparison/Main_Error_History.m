% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 

clc;clear all;close all;
load 'MOTES_Solution.mat'
load 'Analytical_Solution.mat'
xx1 = squeeze(History.Nhist(1,1,1:end));
xx2 = squeeze(History.Nhist(1,2,1:end));
xx3 = squeeze(History.Nhist(1,3,1:end));
% xx4 = squeeze(History.Nhist(1,4,1:end));

yy1 = squeeze(History.Nhist(2,1,1:end));
yy2 = squeeze(History.Nhist(2,2,1:end));
yy3 = squeeze(History.Nhist(2,3,1:end));


%%
figure('Color', 'white')
subplot(2,1,1)
plot(T, xx1(:), T, xx2(:),T, xx3(:),'LineWidth',2)
xlabel('Time','Interpreter','latex')
ylabel('$x$','Interpreter','latex')
set(gca,'fontsize', 15,'linewidth',1.15)
set(gca,'ticklength',1.2*get(gca,'ticklength'))
subplot(2,1,2)
plot(T, yy1(:) , T, yy2(:) ,T,yy3(:),'LineWidth',2)
xlabel('Time','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
legend('Node 1','Node 2','Node 3','Interpreter','latex')
set(gca,'fontsize', 15,'linewidth',1.15)
set(gca,'ticklength',1.2*get(gca,'ticklength'))


figure('Color', 'white')
subplot(2,1,1)
plot(T, zeros(size(x(:,1))), T,  x(:,1),T,x(:,2),'LineWidth',2)
xlabel('Time','Interpreter','latex')
ylabel('$x$','Interpreter','latex')
set(gca,'fontsize', 15,'linewidth',1.15)
set(gca,'ticklength',1.2*get(gca,'ticklength'))
subplot(2,1,2)
plot(T, zeros(size(y(:,1))), T,  y(:,1),T,y(:,2),'LineWidth',2)
xlabel('Time','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
legend('Node 1','Node 2','Node 3','Interpreter','latex')
set(gca,'fontsize', 15,'linewidth',1.15)
set(gca,'ticklength',1.2*get(gca,'ticklength'))
% 

%%
figure('Color', 'white')
subplot(2,1,1)
plot(T, xx1(:)- zeros(size(x(:,1))), T, xx2(:)- x(:,1),T,xx3(:)-x(:,2),'LineWidth',2)
xlabel('Time','Interpreter','latex')
ylabel('$x$','Interpreter','latex')
set(gca,'FontName','Times New Roman','fontsize', 22,'linewidth',1.2)
set(gca,'XMinorTick','on','YMinorTick','on')
set(gca,'ticklength',1.8*get(gca,'ticklength'))
% legend('Node 1','Node 2, 3','Node 4','Interpreter','latex')
legend('Node 1','Node 2','Node 3','Interpreter','latex')
title('Errors in coordinates')
subplot(2,1,2)
plot(T, yy1(:) - zeros(size(y(:,1))), T, yy2(:) - y(:,1),T,yy3(:)-y(:,2),'LineWidth',2)
xlabel('Time','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
set(gca,'FontName','Times New Roman','fontsize', 22,'linewidth',1.2)
set(gca,'XMinorTick','on','YMinorTick','on')
set(gca,'ticklength',1.8*get(gca,'ticklength'))
% legend('Node 1','Node 2, 3','Node 4','Interpreter','latex')
legend('Node 1','Node 2','Node 3','Interpreter','latex')



figure
plot(t,History.bar_len(1,:)-1, t,History.bar_len(2,:)-1,'LineWidth',1)
legend('Bar 1','Bar 2','Interpreter','latex')
xlabel('Time','Interpreter','latex')
ylabel('Error of Bar Length','Interpreter','latex')
set(gca,'fontsize', 15,'linewidth',1.15)
set(gca,'ticklength',1.2*get(gca,'ticklength'))

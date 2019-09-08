function Y = ode4Ck(odefun,tspan,y0,varargin)

% This function is modified from the ode4.m by Kelly Kearney: https://github.com/kakearney/ecosystem-pkg/blob/master/odefixed/ode4.m. 

% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 

%ODE4  Solve differential equations with a non-adaptive method of order 4.
%   Y = ODE4(ODEFUN,TSPAN,Y0) with TSPAN = [T1, T2, T3, ... TN] integrates 
%   the system of differential equations y' = f(t,y) by stepping from T0 to 
%   T1 to TN. Function ODEFUN(T,Y) must return f(t,y) in a column vector.
%   The vector Y0 is the initial conditions at T0. Each row in the solution 
%   array Y corresponds to a time specified in TSPAN.
%
%   Y = ODE4(ODEFUN,TSPAN,Y0,P1,P2...) passes the additional parameters 
%   P1,P2... to the derivative function as ODEFUN(T,Y,P1,P2...). 
%
%   This is a non-adaptive solver. The step sequence is determined by TSPAN
%   but the derivative function ODEFUN is evaluated multiple times per step.
%   The solver implements the classical Runge-Kutta method of order 4.   
%
%   Example 
%         tspan = 0:0.1:20;
%         y = ode4(@vdp1,tspan,[2 0]);  
%         plot(tspan,y(:,1));
%     solves the system y' = vdp1(t,y) with a constant step size of 0.1, 
%     and plots the first component of the solution.   
%

if ~isnumeric(tspan)
  error('TSPAN should be a vector of integration steps.');
end

if ~isnumeric(y0)
  error('Y0 should be a vector of initial conditions.');
end

h = diff(tspan);
if any(sign(h(1))*h <= 0)
  error('Entries of TSPAN are not in order.') 
end  

try
  f0 = feval(odefun,tspan(1),y0,varargin{:});
catch
  msg = ['Unable to evaluate the ODEFUN at t0,y0. ',lasterr];
  error(msg);  
end  

y0 = y0(:);   % Make a column vector.
if ~isequal(size(y0),size(f0))
  error('Inconsistent sizes of Y0 and f(t0,y0).');
end  


varargin=varargin{1,1};
C_b = varargin.tenseg_struct.C_b; len_hat = varargin.tenseg_struct.len_hat;
P = varargin.tenseg_struct.P; D = varargin.tenseg_struct.D;

%     For Reduced Order
[U,SIGMA,V] = svd(P);
rank_SIGMA = rank(SIGMA);
U2 = U(:,rank_SIGMA+1:end);
V1 = V(:,1:rank_SIGMA);
Sigma1 = SIGMA(1:rank_SIGMA,1:rank_SIGMA); 	
X1 = D*V1*Sigma1^-1;

neq = length(y0);
N = length(tspan);
Y = zeros(neq,N);
F = zeros(neq,4);

Y(:,1) = y0;
for i = 2:N
  ti = tspan(i-1);
  hi = h(i-1);
  yi = Y(:,i-1);
  F(:,1) = feval(odefun,ti,yi,varargin);
  F(:,2) = feval(odefun,ti+0.5*hi,yi+0.5*hi*F(:,1),varargin);
  F(:,3) = feval(odefun,ti+0.5*hi,yi+0.5*hi*F(:,2),varargin);  
  F(:,4) = feval(odefun,tspan(i),yi+hi*F(:,3),varargin);
  Y(:,i) = yi + (hi/6)*(F(:,1) + 2*F(:,2) + 2*F(:,3) + F(:,4));
%%
X2 = Y(1:numel(y0)/2,i);
X2d = Y(numel(y0)/2+1:end,i);

X2 = reshape(X2,3,size(U2,2));
X2d = reshape(X2d,3,size(U2,2));

N = [X1 X2]*U';
Nd = [0*X1 X2d]*U';

%%%%bar length correction
[N,Nd]=bar_length_correction_ClassK(N,Nd,C_b,P,D,len_hat);
X2 = N*U2;
X2d = Nd*U2;
Y(:,i)=[X2(:);X2d(:)];
end
Y = Y.';

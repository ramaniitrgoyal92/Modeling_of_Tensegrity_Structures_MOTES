function [Force_den,Support_Load] = tenseg_equilibrium(tenseg_struct)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 
%
% [Force_den,Support_Load] = tenseg_equilibrium(tenseg_struct) Solves the 
% statics equation NK=W
%
%
% Inputs:  
%	tenseg_struct: tensegrity structure definition data structure
%		containing bar and string connectivity matrices
% Outputs:
%	Force_den: the Force_den in the strings and bars
%	Support_Load: the Support_Load at the constriant nodes

% % Load in structure definitions
N = tenseg_struct.N;
C_s = tenseg_struct.C_s;
C_b = tenseg_struct.C_b;
W = tenseg_struct.W;
Pinned_nodes = tenseg_struct.Pinned_nodes;

% clear all; clc;
% [N,C_b,C_s] = tenseg_prismplate(1,1);
% 
% W = N;
% Pinned_nodes = [2 1 0 1;4 0 0 1];

if ~isempty(C_b)
    for i = 1:size(N,2)
    E = [C_s'*diag(C_s(:,i)) -C_b'*diag(C_b(:,i))];
    K(3*i-2:3*i,:) = N*E;
    end
else
    for i = 1:size(N,2)
    E = C_s'*diag(C_s(:,i));
    K(3*i-2:3*i,:) = N*E;
    end
end

K1 = K;
W_vec = reshape(W,[],1);

if ~isempty(Pinned_nodes)
% Pinned Nodes
remove = [];
j = 1;
nc = Pinned_nodes(:,1)';
for i = 1:size(nc,2)
    if Pinned_nodes(i,2) == 1
        K2(j,:) = K(3*nc(i)-2,:);
        remove = [remove 3*nc(i)-2];
        j = j+1;
    end
    if Pinned_nodes(i,3) == 1
        K2(j,:) = K(3*nc(i)-1,:);
        remove = [remove 3*nc(i)-1];
        j = j+1;
    end
    if Pinned_nodes(i,4) == 1
        K2(j,:) = K(3*nc(i),:);
        remove = [remove 3*nc(i)];
        j = j+1;
    end
end


% Defined Loads
K1(remove,:) = []; 
W_vec(remove,:) = [];

end

% size(K1)
% rank(K1)
% Force_den = null(K1)
% null(K1);
% rank([K1 W_vec])

%% -------------- Steel Properties for Bar -----------
% E_b = 200e09;
% rho_b = 8000;
% sigma_b = 300e06;

% -------------- Carbon fibre Properties for Bar -----------
E_b = 230e09;
rho_b = 1400;
sigma_b = 3500e06;

%---------------- UHMWPE Properties for Bar -------------------
% E_b = 120e09; 
% rho_b = 970;
% sigma_b = 2.7e09;

%---------------- Aluminum for Bar -------------------
% E_b = 60e09;
% rho_b = 2700;
% sigma_b = 110e06;


%% -------------- Steel Properties for Strings -----------
% E_s = 200e09;
% rho_s = 8000;
% sigma_s = 300e06;

%------------------ UHMWPE Properties for Strings -----------
E_s = 120e09; 
rho_s = 970;
sigma_s = 2.7e09;


%---------------- Aluminum for Strings -------------------
% E_s = 60e09;
% rho_s = 2700;
% sigma_s = 110e06;

%% Final Optimization

B = N*C_b';
S = N*C_s';

Bar_len = diag(sqrt(B'*B));
String_len = diag(sqrt(S'*S));

mass_gamma = rho_s/sigma_s*String_len.*String_len;
mass_lambda = rho_b/sigma_b*Bar_len.*Bar_len;


options = optimoptions('linprog','Algorithm','interior-point-legacy','maxiterations',2000);
% options = optimoptions('linprog','Algorithm','dual-simplex');%,'maxiterations',2000);
% Force_den = linprog(ones(1,size(C_s,1)+size(C_b,1)),[],[],K1,W_vec,zeros(1,size(C_s,1)),[],options);
% Force_den = linprog(ones(1,size(C_s,1)+size(C_b,1)),[],[],K1,W_vec,5e4*ones(1,size(C_s,1)),[],options);
% Force_den = linprog([mass_gamma;mass_lambda],[],[],K1,W_vec,1e1*ones(1,size(C_s,1)),[],options);
% Force_den = linprog([],[],[],K1,W_vec,zeros(1,size(C_s,1)),[],options);
% Force_den = linprog([],[],[],K1,zeros(size(W_vec)),.01*ones(1,size(C_s,1)),[],options);
% Force_den = linprog([],[],[],K1,W_vec,.01*ones(1,size(C_s,1)),[],options);
Force_den = linprog([mass_gamma;mass_lambda],[],[],K1,W_vec,zeros(1,size(C_s,1)),[],options);
% Force_den = pinv(K1)*W_vec;


Mass_string = mass_gamma'*Force_den(1:size(C_s,1),1);
Mass_bar = mass_lambda'*Force_den(1+size(C_s,1):end,1);

Mass_bar_buckling = ((2/sqrt(pi*E_b))*rho_b*Bar_len.*Bar_len)'*sqrt(Force_den(1+size(C_s,1):end,1).*Bar_len);

% K1*Force_den
if isempty(Pinned_nodes)
    Support_Load = 0;
else
    Support_Load = K2*Force_den;
end

function [Force_den,Support_Load,LABEL,MIN_MASS,Mass_bar,Mass_string,BarWidth_in,StringWidth_in,Loop] = tenseg_equilibrium_minimal_mass(tenseg_struct)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 


% INPUT: STRUCTURES
% OUTPUT1: FORCE_DENSITY;
% OUTPUT2: SUPPORT_LOAD;
% OUTPUT3: LABEL; IF THE ENTRIES ARE 1, IT'S BUCKLING; OTHERWISE IT'S
% YIELDING

% Load in structure definitions
N = tenseg_struct.N;
C_s = tenseg_struct.C_s;
C_b = tenseg_struct.C_b;
W = tenseg_struct.W;
% bar_failure=tenseg_struct.bar_failure; % 'yielding' or 'yielding_n_buckling'
bar_material=tenseg_struct.bar_material;
string_material=tenseg_struct.string_material;
Pinned_nodes = tenseg_struct.Pinned_nodes;

for i = 1:size(N,2)
    E = [C_s'*diag(C_s(:,i)) -C_b'*diag(C_b(:,i))];
    K(3*i-2:3*i,:) = N*E;
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

size(K1);
rank(K1);
% null(K1)
% null(K1')
rank([K1 W_vec]);


switch bar_material
    case 'Steel'
        % -------------- Steel Properties for Bar -----------
        E_b = 200e09;
        rho_b = 8000;
        sigma_b = 300e06;
    case 'UHMWPE'
        %---------------- UHMWPE Properties for Bar ----------
        E_b = 120e09;
        rho_b = 970;
        sigma_b = 2.7e09;
    case 'Aluminum'
        %---------------- Aluminum for Bar -------------------
        E_b = 60e09;
        rho_b = 2700;
        sigma_b = 110e06;
    otherwise
        disp('Edit the material database')
end

switch string_material
    case 'Steel'
        %% -------------- Steel Properties for Strings -----------
        E_s = 200e09;
        rho_s = 8000;
        sigma_s = 300e06;
    case 'UHMWPE'
        %------------------ UHMWPE Properties for Strings ---------
        E_s = 120e09;
        rho_s = 970;
        sigma_s = 2.7e09;
    case 'Aluminum'
        %---------------- Aluminum for Strings -------------------
        E_s = 60e09;
        rho_s = 2700;
        sigma_s = 110e06;
    otherwise
        disp('Edit the material database')
end
%% Final Optimization
B = N*C_b';
S = N*C_s';
ns=size(S,2);
nb=size(B,2);

% Square of bar, stirng length
Bar_len = diag(B'*B);
String_len = diag(S'*S);

mass_gamma = rho_s/sigma_s*String_len;
mass_lambda = rho_b/sigma_b*Bar_len;
mass_member=[mass_gamma;mass_lambda];

mass_gamma_buckle = rho_s/sigma_s*String_len;
mass_lambda_buckle = 2*rho_b/sqrt(pi*E_b)*Bar_len.^(5/4);
mass_member_buckle=[mass_gamma_buckle;mass_lambda_buckle];

%% Assume that all bars will buckle, this leads to all labels for bar to be
% 1, otherwise it is zero.
Force_den=cell(50,1);

LABEL=ones(nb,1);
x0=1*10^-2.*ones(ns+nb,1);
options = optimoptions('fmincon','Algorithm','sqp'); 
Force_den{1}= fmincon(@(x) nonlinear_force_density_minimal_mass(x,ns,LABEL,...
    mass_member_buckle,mass_member),x0,[],[],K1,W_vec,zeros(1,size(C_s,1)))

% Loop=1;

%% Check if the bar force densities satisfy the buckling assumption
lambda_test=Force_den{1}(length(S(1,:))+1:end);
for i=1:nb
    if lambda_test(i)>=4*sigma_b^2*sqrt(Bar_len(i))/(pi*E_b)
        LABEL(i)=0; % MEANS THAT ITS YIELDING INSTEAD OF BUCKLING
    else
        LABEL(i)=LABEL(i);  % MEANS THAT ITS BUCKLING AS WE ASSUME
    end
end

Loop=2;
x0=1*10^-2.*ones(ns+nb,1);
options = optimoptions('fmincon','Algorithm','sqp');
Force_den{Loop}= fmincon(@(x) nonlinear_force_density_minimal_mass(x,ns,LABEL,...
    mass_member_buckle,mass_member),x0,[],[],K1,W_vec,zeros(1,size(C_s,1)))

%% While loop for converging

while abs(Force_den{Loop}-Force_den{Loop-1})>10^-2
    
    lambda_test=Force_den{Loop}(length(S(1,:))+1:end);
    for i=1:nb
        if lambda_test(i)>=4*sigma_b^2*sqrt(Bar_len(i))/(pi*E_b)
            LABEL(i)=0; % MEANS THAT ITS YIELDING INSTEAD OF BUCKLING
        else
            LABEL(i)=LABEL(i);  % MEANS THAT ITS BUCKLING AS WE ASSUME
        end
    end
    
    Loop=Loop+1;
    x0=1*10^-2.*ones(ns+nb,1);
    options = optimoptions('fmincon','Algorithm','sqp')
    Force_den{Loop}= fmincon(@(x) nonlinear_force_density_minimal_mass(x,ns,LABEL,...
        mass_member_buckle,mass_member),x0,[],[],K1,W_vec,zeros(1,size(C_s,1)))
    
end

%% Output the mass

MM=[];Mass_bar=[];Mass_string=[];


MM(1:ns)=mass_member_buckle(1:ns).*Force_den{Loop}(1:ns); %%%%%%%%%%%%%%%%%%%%%%%%%%

for i=ns+1:length(mass_member_buckle)
    if LABEL(i-ns)==1;
        MM(i)=mass_member_buckle(i)*sqrt(Force_den{Loop}(i));
    else
        MM(i)=mass_member(i)*Force_den{Loop}(i);
    end
end

Mass_string = MM(1:ns);%%%%%%%%%%%%%%%%%%%%%
Mass_bar = MM(ns+1:end);
[MIN_MASS,FF] = nonlinear_force_density_minimal_mass(Force_den{Loop},ns,LABEL,mass_member_buckle,mass_member);
%% Here, we wants to show in plots that bars ans strings thinkness is proportional to their actural radius.

r_bar=[];r_string=[];BarWidth_in=[];StringWidth_in=[];
% Radius of Bars 
r_bar=(Mass_bar'./(rho_b*pi*Bar_len.^(1/2))).^(1/2);
% Radius of Strings
r_string=(Mass_string'./(rho_s*pi*String_len.^(1/2))).^(1/2);

if max(r_bar)-min(r_bar) < 10^-6 % If bars with same radius, set it as default line width 8.
	BarWidth_in = 8*ones(size(r_bar));
else
    % Bars with least radius, line width is 2, Bars with largest radius,
    % line width is 8.
    BarWidth_in=2+(8-2)*(r_bar-min(r_bar))/(max(r_bar)-min(r_bar)); 
end 

% q=1, k=4, r_string too small, so set to mim
if max(r_string) - min(r_string) < 10^-6
    StringWidth_in=2*ones(size(r_string)); % If strings with almost same radius, set it as default line width 2.
else
    % Strings with least radius, line width is 2, Strings with largest radius,
    % line width is 8.
    StringWidth_in=2+(8-2)*(r_string-min(r_string))/(max(r_string)-min(r_string));
end
%%
if isempty(Pinned_nodes)
    Support_Load = 0;
else
    Support_Load = K2*Force_den{Loop};
end

end

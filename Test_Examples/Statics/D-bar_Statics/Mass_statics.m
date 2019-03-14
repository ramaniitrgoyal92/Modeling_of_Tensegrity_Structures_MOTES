%% Total mass calculator
function [total_mass,bar_mass,string_mass]=Mass_statics(tenseg_struct,Force_den)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 

% Load in structure definitions
N = tenseg_struct.N;
C_s = tenseg_struct.C_s;
C_b = tenseg_struct.C_b;
S = N*C_s'; B = N*C_b';
bar_failure=tenseg_struct.bar_failure; % 'yielding' or 'yielding_n_buckling'
bar_material=tenseg_struct.bar_material;
string_material=tenseg_struct.string_material;
gamma=Force_den(1:length(S(1,:))); % take gamma out of the force density vector 
lambda=Force_den(length(S(1,:))+1:end); % take lambda out of the force density vector 
%% get the right material for the structure system
% bar material
switch bar_material
    case 'Steel'
        %-------------- Steel Properties for Bar -------------
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
% string material
switch string_material
    case 'Steel'
        %-------------- Steel Properties for Strings -----------
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
%% Solve for structure mass by force densities subject to yielding
string_mass = rho_s/sigma_s*trace(S*diag(gamma)*S'); % string yielding
% consider two failure modes cases for bars
save text.mat
switch bar_failure
    case 'yielding'
        bar_mass_yld = rho_b/sigma_b*trace(B*diag(lambda)*B'); % bar yielding
        bar_mass = bar_mass_yld;
    case 'yielding_n_buckling'
        bar_mass_max = 0;
        for j = 1 : length(B(1,:)) % number of bars
            bar_mass_yld = rho_b/sigma_b*norm(B(:,j))^2*lambda(j); % bar yielding
            bar_mass_bkl = 2*rho_b/sqrt(pi*E_b)*sqrt(lambda(j))*power(norm(B(:,j)),2.5); % bar buckling
            bar_mass_max = bar_mass_max + max(bar_mass_yld,bar_mass_bkl); % take the bigger mass
            if bar_mass_bkl > bar_mass_yld
                disp(strcat('Buckling dominant: Bar',32,num2str(j))) % 32 is a space note
            else
                disp(strcat('Yielding dominant: Bar',32,num2str(j))) % 32 is a space note
            end
        end
        bar_mass = bar_mass_max;
end
total_mass = string_mass + bar_mass;

end
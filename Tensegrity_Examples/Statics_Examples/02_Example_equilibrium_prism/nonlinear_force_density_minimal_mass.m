function [F,FF]=nonlinear_force_density_minimal_mass(x,ns,LABEL,mass_member_buckle,mass_member)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 

% This function is calculation mass of each structure member and the total
% mass.

% Inputs:
%	x: force densities in strings and bars
%   ns: number of strings
%   LABEL: a vector denotes bar is yielding or buckling
%   mass_member_buckle: material constants for buckling calculation
%   mass_member: material constants for yielding calculation
% Outputs:
%	F: mass of each structure element
%	FF: Total mass of all the bars and strings

FF=[];
FF(1:ns)=mass_member_buckle(1:ns).*x(1:ns); % Mass of strings subject to yileding

% Mass of bars, the calculation is based on the labels (denoting bar is yielding or buckling)
for i=ns+1:length(mass_member_buckle)
    if LABEL(i-ns)==1
        FF(i)=mass_member_buckle(i)*sqrt(x(i)); % Mass subject to Buckling
    else
        FF(i)=mass_member(i)*x(i); % Mass subject to Yielding
    end
end

F=sum(FF); % Total Mass
end
 
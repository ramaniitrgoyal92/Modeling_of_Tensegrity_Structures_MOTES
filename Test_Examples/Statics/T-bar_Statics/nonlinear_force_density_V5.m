
function [F,FF]=nonlinear_force_density_V5(x,ns,LABEL,mass_member_new,mass_member)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 
FF=[];
FF(1:ns)=mass_member_new(1:ns).*x(1:ns);

for i=ns+1:length(mass_member_new)
    if LABEL(i-ns)==1
        FF(i)=mass_member_new(i)*sqrt(x(i));
    else
        FF(i)=mass_member(i)*x(i);
    end
end

F=sum(FF);
end
 
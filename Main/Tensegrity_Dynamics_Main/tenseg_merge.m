function [N,CST,CBT] = tenseg_merge(N1,cst1,cbt1,N2,cst2,cbt2)

% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 

% [N,CST,CBT] = tenseg_merge(N1,cst1,cbt1,N2,cst2,cbt2) Merges two tensegrity 
% structures into one.

%
% Inputs:
%	N1: node matrix (3 x n1 array for n1 nodes) of structure 1
%	cst1: string connectivity matrix (alpha1 x n1 array for alpha1 strings) of structure 1
%	cbt1: bar connectivity matrix (beta1 x n1 array for beta1 bars) of structure 1
%	N2: node matrix (3 x n2 array for n2 nodes) of structure 2
%	cst2: string connectivity matrix (alpha2 x n2 array for alpha2 strings) of structure 2
%	cbt2: bar connectivity matrix (beta2 x n12 array for beta12 bars) of structure 2
%
% Outputs:
%	N: merged node matrix
%	CST: merged bar connectivity matrix
%	CBT: merged string connectivity matrix





pos1 = [];
pos2 = [];
for i = 1:size(N1,2)
    for j = 1:size(N2,2)
        if norm(N1(:,i)-N2(:,j)) <= 1e-4
            pos1 = [pos1 i];
            pos2 = [pos2 j];
        end
    end
end
N2(:,pos2) = [];
N = [N1 N2];


% ==================== String Connectivity
cst1 = [cst1;zeros(size(N2,2),size(cst1,2))];
cst2_temp = cst2(pos2,:);
cst2(pos2,:) = [];
cst2_upper = zeros(size(N1,2),size(cst2,2));

for i = 1:size(pos1,2)
    cst2_upper(pos1(1,i),:) = cst2_temp(i,:);
end


cst2 = [cst2_upper;cst2];

CST = [cst1 cst2];

pos3 = [];
for i = 1:size(CST,2)
    for j = i+1:size(CST,2)
        if norm(CST(:,i)-CST(:,j)) == 0 || norm(CST(:,i)-CST(:,j)) == sqrt(8)
            pos3 = [pos3 i];
        end
    end
end

CST(:,pos3)=[];

        

% ==================== Bar Connectivity

cbt1 = [cbt1;zeros(size(N2,2),size(cbt1,2))];
cbt2_temp = cbt2(pos2,:);
cbt2(pos2,:) = [];
cbt2_upper = zeros(size(N1,2),size(cbt2,2));

for i = 1:size(pos1,2)
    cbt2_upper(pos1(1,i),:) = cbt2_temp(i,:);
end


cbt2 = [cbt2_upper;cbt2];

CBT = [cbt1 cbt2];

pos4 = [];
for i = 1:size(CBT,2)
    for j = i+1:size(CBT,2)
        if norm(CBT(:,i)-CBT(:,j)) == 0 || norm(CBT(:,i)-CBT(:,j)) == sqrt(8)
            pos4 = [pos4 i];
        end
    end
end

CBT(:,pos4)=[];
end
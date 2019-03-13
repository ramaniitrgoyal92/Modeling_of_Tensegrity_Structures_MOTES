function [N_Cor,Ndot_Cor]=bar_length_correction_Class1(N,Ndot,C_b,bar_len_const_hat)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 

% [N_Cor,Ndot_Cor]=bar_length_correction_Class1(N,Ndot,C_b,bar_len_const_hat)
% is to correct the bar length errors due to numerical simulations
% errors for Class1 structures.

%
% Inputs:
%	N: current nodel matrix
%	Ndot: current N_dot values
%	C_b: bar connectivity matrix
%   bar_len_const_hat: constant bar length matrix

% Outputs:
%	N_Cor: corrected nodel matrix
%	Ndot_Cor: corrected N_dot matrix



C_r = 1/2*abs(C_b);
q=.50;

B=N*C_b';Bdot=Ndot*C_b';
R=N*C_r';Rdot=Ndot*C_r';

for j=1:size(B,2)
    
    a0=(q*bar_len_const_hat(j,j))^2*norm(Bdot(:,j))^2*((B(:,j)'*Bdot(:,j))^2-norm(B(:,j))^2*norm(Bdot(:,j))^2);
    a1=2*(q*bar_len_const_hat(j,j))^2*((B(:,j)'*Bdot(:,j))^2-norm(B(:,j))^2*norm(Bdot(:,j))^2);
    a2=norm(Bdot(:,j))^4-(q*bar_len_const_hat(j,j))^2*norm(B(:,j))^2;
    a3=2*norm(Bdot(:,j))^2;
    x=roots([1 a3 a2 a1 a0]);
    x= x(imag(x)==0);
    J=zeros(1,size(x,1));
    for i=1:size(x,1)
        v=q*bar_len_const_hat(j,j)*((x(i)*eye(3)+Bdot(:,j)*Bdot(:,j)')\B(:,j));
        p=bar_len_const_hat(j,j)*v-B(:,j);
        r=-v*v'*Bdot(:,j);
        J(i)=q*norm(p)^2+norm(r)^2;
    end
    [~,min_index]=min(J);
    v=q*bar_len_const_hat(j,j)*((x(min_index)*eye(3)+Bdot(:,j)*Bdot(:,j)')\B(:,j));
    p=bar_len_const_hat(j,j)*v-B(:,j);
    r=-v*v'*Bdot(:,j);
    B(:,j)=B(:,j)+p;
    Bdot(:,j)=Bdot(:,j)+r;
end

    Q=[B R];
    Qdot=[Bdot Rdot];
    N_Cor=Q/[C_b' C_r'];
    Ndot_Cor=Qdot/[C_b' C_r'];


if norm(C_b(:,size(C_b,2)))==0
    for k=1:size(C_b,2)-2*size(B,2)
        N_Cor(:,2*size(B,2)+k)=N(:,2*size(B,2)+k);
        Ndot_Cor(:,2*size(B,2)+k)=Ndot(:,2*size(B,2)+k);
    end
end
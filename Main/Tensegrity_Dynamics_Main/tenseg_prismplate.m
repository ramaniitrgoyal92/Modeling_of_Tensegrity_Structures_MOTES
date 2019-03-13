function [N,C_b,C_s] = tenseg_prismplate(q,Height,Dia)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 

% [N,C_b,C_s] = TENSEG_PRISMPLATE(q,Height) generates node and connectivity
% matrices for a prism-based tensegrity plate with variable complexity and
% height.
%
% Inputs:
%	q: complexity number
%	Height: height of plate
% 
% Outputs:
%	N: node matrix
%	C_b: bar connectivity matrix
%	C_s: string connectivity matrix

%***************************************************************************
% 1.form the connectivity matrix and the initial node location and choose
% the initial preforce.
%*******************************************************************************
% 1.1 input parameters
% Dia=1;%cover diameter
Area=3*(Dia/2)^2*cos(pi/6);%cover area
% Height=0.01;%the height of plate
% q=5;%the complexity q

% 1.2. connectivity matrix
[CBT,CST,CRT]=connectivitymatrix(q);
% 1.3 initial node location Node_I and bar string matrix
[Node_I,r0]=nodematrix(q,Height,Area);%node matrix
nb=size(CBT,2);ns=size(CST,2);
%plot
[rows,cols]=find(CST);[rowb,colb]=find(CBT);
% figure;
% for j=1:nb
% plot3(Node_I(1,rowb(2*j-1:2*j)),Node_I(2,rowb(2*j-1:2*j)),Node_I(3,rowb(2*j-1:2*j)),'linewidth',2,'color','k'); hold on;
% end
% for j=1:ns
% plot3(Node_I(1,rows(2*j-1:2*j)),Node_I(2,rows(2*j-1:2*j)),Node_I(3,rows(2*j-1:2*j)),'linewidth',1,'color','r'); hold on;
% end
% xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
% grid on;
% axis equal;
% drawnow;
% Node_I
% CBT
% CST

N = Node_I;
C_s = CST';
C_b = CBT';

end

%*****************************
% sub function
%*********************************
%*******************************
% 1. connectivity matrix
%********************************
function [CBT,CST,CRT]=connectivitymatrix(q) 
nu=1+3*q*(q-1);%number of units
nn=6*nu;%number of nodes
ns=(1+3*(q-1)*(q-2))*15+(6*(q-1)-3)*13+3*11;%number of strings
% 1) for bars
I3=speye(3);Ib=[-I3;I3];Ir=0.5*[I3;I3];In=speye(nu);
CBT=kron(In,Ib);CRT=kron(In,Ir);
         
% 2) for strings
e0=sparse(zeros(6,1));ee=sparse(eye(6));e1=ee(:,1);e2=ee(:,2);e3=ee(:,3);e4=ee(:,4);e5=ee(:,5);e6=ee(:,6);
% for s(0,0,0)
Q00=[-e1 e2 -e2 e3 -e3 e1 -e4 e5 -e5 e6 -e6 e4 e3-e4 e1-e5 e2-e6];
Q01=[e3 -e3 e0 e0 e0 e0 e0 e0 e4 -e4 e0 e0 e0 e0 e0];
Q02=[e0 e0 e3 -e3 e0 e0 e0 e0 e0 e0 e4 -e4 e0 e0 e0];
Q03=[e0 e0 e0 e0 e3 -e3 e4 -e4 e0 e0 e0 e0 e0 e0 e0];
% for s(j,0,0)
Qj0_tl=Q00;Qj1_tl=Q01;Qj3_tl=Q03;
Q12_tl=[e0 e0 e1 -e1 e0 e0 e0 e0 e0 e0 e5 -e5 e0 e0 e0];
Q22_tl=[e0 e0 e2 -e2 e0 e0 e0 e0 e0 e0 e6 -e6 e0 e0 e0];
Q32_tl=[e0 e0 e3 -e3 e0 e0 e0 e0 e0 e0 e4 -e4 e0 e0 e0];
% for s(j,0,k)
Qj0_t=Q00;Qj1_t=Q01;Qj2_t=Q02;
Qj3_t=[e0 e0 e0 e0 e2 -e2 e6 -e6 e0 e0 e0 e0 e0 e0 e0];
% for s(j,0,q-1)
Qj0_tr=[Q00(:,[3:8 11:15]) e2-e1 e6-e5];Qj2_tr=[Q02(:,[3:8 11:15]) e0 e0];Qj3_tr=[Qj3_t(:,[3:8 11:15]) e0 e0];
% for s(j,i,0)
Qj0_ml=Q00;Qj1_ml=Q01;Qj3_ml=Q03;
Qj2_ml=[e0 e0 e1 -e1 e0 e0 e0 e0 e0 e0 e5 -e5 e0 e0 e0];
% for s(j,i,k)
Qj0_m=Q00;Qj1_m=Q01;Qj2_m=Qj2_ml;Qj3_m=Qj3_t;
% for s(j,i,q-1)
Qj0_mr=Qj0_tr;Qj2_mr=[Qj2_ml(:,[3:8 11:15]) e0 e0];Qj3_mr=Qj3_tr;
% for s(j,q-2,0)
Qj0_bl=[Q00(:,[3 4 11:15]) e2-e1 e1-e3 e5-e4 e6-e5];Qj2_bl=[Qj2_ml(:,[3 4 11:15]) e0 e0 e0 e0];
% for s(j,q-2,k)
Qj0_b=Qj0_tr;Qj2_b=Qj2_mr;Qj3_b=Qj3_tr;
% for s(j,q-2,q-1)
Qj0_br=Qj0_tr;Qj2_br=Qj2_mr;Qj3_br=Qj3_tr;


CST=sparse(zeros(nn,ns)); 
if q==1
    CST=[e3-e4 e1-e5 e2-e6 e2-e1 e3-e2 e1-e3 e5-e4 e6-e5 e4-e6];
    count_col=9;
    Cou_col=[1;9];
else
count_col=15*ones(1,nu);%count_col is a vector store the string number of each units,which is the column number of CST
count_row=6*ones(1,nu);%count_row is a vector store the node number of each units,which is the row number of CST
for j=1:3
    for i=0:q-2
        for k=0:q-1
            nou=2+(q-1)*q*(j-1)+q*i+k;
            if (i~=q-2)&&(k==q-1)
                count_col(nou)=13;
            elseif (i==q-2)&&(k==0)
                count_col(nou)=11;
            elseif (i==q-2)&&(k~=0)
                count_col(nou)=13;
            end         
        end
    end
end
Cou_col=zeros(2,nu);%the column block range of CST corresponding to each unit 
Cou_row=zeros(2,nu);%the row block range of CST corresponding to each unit
for n=1:nu
   Cou_col(1,n)=sum(count_col(1,1:n))-count_col(1,n)+1; 
   Cou_col(2,n)=sum(count_col(1,1:n));
   Cou_row(1,n)=sum(count_row(1,1:n))-count_row(1,n)+1; 
   Cou_row(2,n)=sum(count_row(1,1:n));
end
% FOR(0,0,0)
nc=1;nr0=nc;nr1=3;nr2=3+(q-1)*q;nr3=3+(q-1)*q*2;
CST(Cou_row(1,nr0):Cou_row(2,nr0),Cou_col(1,nc):Cou_col(2,nc))=Q00;
CST(Cou_row(1,nr1):Cou_row(2,nr1),Cou_col(1,nc):Cou_col(2,nc))=Q01;
CST(Cou_row(1,nr2):Cou_row(2,nr2),Cou_col(1,nc):Cou_col(2,nc))=Q02;
CST(Cou_row(1,nr3):Cou_row(2,nr3),Cou_col(1,nc):Cou_col(2,nc))=Q03;
if q==2
    for j=1:3
    %for(j,0,0)
    if j==1
        Qj2_tl=Q12_tl;
    elseif j==2
        Qj2_tl=Q22_tl;
    elseif j==3
        Qj2_tl=Q32_tl;
    end
    nc=(j-1)*(q-1)*q+2;nr0=nc;nr2=1;
    CST(Cou_row(1,nr0):Cou_row(2,nr0),Cou_col(1,nc):Cou_col(2,nc))=[Q00(:,[3 4 11:15]) e2-e1 e1-e3 e5-e4 e6-e5];
    CST(Cou_row(1,nr2):Cou_row(2,nr2),Cou_col(1,nc):Cou_col(2,nc))=[Qj2_tl(:,[3 4 11:15]) e0 e0 e0 e0];
    %for(j,0,1)
    nc=(j-1)*(q-1)*q+q+1;nr0=nc;nr3=nr0-1;
    if j==3
        nr2=(q-2)*q+2;
    else
        nr2=j*(q-1)*q+(q-2)*q+2;
    end
    CST(Cou_row(1,nr0):Cou_row(2,nr0),Cou_col(1,nc):Cou_col(2,nc))=Qj0_tr;
    CST(Cou_row(1,nr2):Cou_row(2,nr2),Cou_col(1,nc):Cou_col(2,nc))=Qj2_tr;
    CST(Cou_row(1,nr3):Cou_row(2,nr3),Cou_col(1,nc):Cou_col(2,nc))=Qj3_tr;
    end

else
for j=1:3
    %FOR(j,0,0)
    if j==1
        Qj2_tl=Q12_tl;nr3=4+(q-1)*q*2;
    elseif j==2
        Qj2_tl=Q22_tl;nr3=4;
    elseif j==3
        Qj2_tl=Q32_tl;nr3=4+(q-1)*q;
    end
    nc=(j-1)*(q-1)*q+2;nr0=nc;nr1=(j-1)*(q-1)*q+3+q;nr2=1;
    CST(Cou_row(1,nr0):Cou_row(2,nr0),Cou_col(1,nc):Cou_col(2,nc))=Qj0_tl;
    CST(Cou_row(1,nr1):Cou_row(2,nr1),Cou_col(1,nc):Cou_col(2,nc))=Qj1_tl;
    CST(Cou_row(1,nr2):Cou_row(2,nr2),Cou_col(1,nc):Cou_col(2,nc))=Qj2_tl;
    CST(Cou_row(1,nr3):Cou_row(2,nr3),Cou_col(1,nc):Cou_col(2,nc))=Qj3_tl;
    %FOR(j,0,k)
    for k=1:q-2
        nc=(j-1)*(q-1)*q+2+k;nr0=nc;nr1=nr0+q+1;nr3=nr0-1;
        if j==3
            nr2=(k-1)*q+2;
        else
            nr2=j*(q-1)*q+(k-1)*q+2;
        end
        CST(Cou_row(1,nr0):Cou_row(2,nr0),Cou_col(1,nc):Cou_col(2,nc))=Qj0_t;
        CST(Cou_row(1,nr1):Cou_row(2,nr1),Cou_col(1,nc):Cou_col(2,nc))=Qj1_t;
        CST(Cou_row(1,nr2):Cou_row(2,nr2),Cou_col(1,nc):Cou_col(2,nc))=Qj2_t;
        CST(Cou_row(1,nr3):Cou_row(2,nr3),Cou_col(1,nc):Cou_col(2,nc))=Qj3_t;
    end
    %FOR(j,0,q-1)
    nc=(j-1)*(q-1)*q+q+1;nr0=nc;nr3=nr0-1;
    if j==3
        nr2=(q-2)*q+2;
    else
        nr2=j*(q-1)*q+(q-2)*q+2;
    end
    CST(Cou_row(1,nr0):Cou_row(2,nr0),Cou_col(1,nc):Cou_col(2,nc))=Qj0_tr;
    CST(Cou_row(1,nr2):Cou_row(2,nr2),Cou_col(1,nc):Cou_col(2,nc))=Qj2_tr;
    CST(Cou_row(1,nr3):Cou_row(2,nr3),Cou_col(1,nc):Cou_col(2,nc))=Qj3_tr;
    
    
    for i=1:q-3
        %for(j,i,0)
        nc=(j-1)*(q-1)*q+i*q+2;nr0=nc;nr1=nr0+q+1;nr2=nr0-q;
        if j==1
            nr3=(q-1)*q*2+i+4;
        else
            nr3=(j-2)*(q-1)*q+i+4;
        end
        CST(Cou_row(1,nr0):Cou_row(2,nr0),Cou_col(1,nc):Cou_col(2,nc))=Qj0_ml;
        CST(Cou_row(1,nr1):Cou_row(2,nr1),Cou_col(1,nc):Cou_col(2,nc))=Qj1_ml;
        CST(Cou_row(1,nr2):Cou_row(2,nr2),Cou_col(1,nc):Cou_col(2,nc))=Qj2_ml;
        CST(Cou_row(1,nr3):Cou_row(2,nr3),Cou_col(1,nc):Cou_col(2,nc))=Qj3_ml;
        %for(j,i,k)
        for k=1:q-2
            nc=(j-1)*(q-1)*q+i*q+k+2;nr0=nc;nr1=nr0+q+1;nr2=nr0-q;nr3=nr0-1;
            CST(Cou_row(1,nr0):Cou_row(2,nr0),Cou_col(1,nc):Cou_col(2,nc))=Qj0_m;
            CST(Cou_row(1,nr1):Cou_row(2,nr1),Cou_col(1,nc):Cou_col(2,nc))=Qj1_m;
            CST(Cou_row(1,nr2):Cou_row(2,nr2),Cou_col(1,nc):Cou_col(2,nc))=Qj2_m;
            CST(Cou_row(1,nr3):Cou_row(2,nr3),Cou_col(1,nc):Cou_col(2,nc))=Qj3_m;            
        end
        %for(j,i,q-1)
        nc=(j-1)*(q-1)*q+i*q+q+1;nr0=nc;nr2=nr0-q;nr3=nr0-1;
        CST(Cou_row(1,nr0):Cou_row(2,nr0),Cou_col(1,nc):Cou_col(2,nc))=Qj0_mr;
        CST(Cou_row(1,nr2):Cou_row(2,nr2),Cou_col(1,nc):Cou_col(2,nc))=Qj2_mr;
        CST(Cou_row(1,nr3):Cou_row(2,nr3),Cou_col(1,nc):Cou_col(2,nc))=Qj3_mr;
    end
    
    %for(j,q-2,0)
    nc=(j-1)*(q-1)*q+(q-2)*q+2;nr0=nc;nr2=nr0-q;
    CST(Cou_row(1,nr0):Cou_row(2,nr0),Cou_col(1,nc):Cou_col(2,nc))=Qj0_bl;
    CST(Cou_row(1,nr2):Cou_row(2,nr2),Cou_col(1,nc):Cou_col(2,nc))=Qj2_bl;
    %for(j,q-2,k)
    for k=1:q-2
        nc=(j-1)*(q-1)*q+(q-2)*q+k+2;nr0=nc;nr2=nr0-q;nr3=nr0-1;
        CST(Cou_row(1,nr0):Cou_row(2,nr0),Cou_col(1,nc):Cou_col(2,nc))=Qj0_b;
        CST(Cou_row(1,nr2):Cou_row(2,nr2),Cou_col(1,nc):Cou_col(2,nc))=Qj2_b;
        CST(Cou_row(1,nr3):Cou_row(2,nr3),Cou_col(1,nc):Cou_col(2,nc))=Qj3_b;        
    end
    %for(j,q-2,q-1)
    nc=(j-1)*(q-1)*q+(q-2)*q+q+1;nr0=nc;nr2=nr0-q;nr3=nr0-1;
    CST(Cou_row(1,nr0):Cou_row(2,nr0),Cou_col(1,nc):Cou_col(2,nc))=Qj0_br;
    CST(Cou_row(1,nr2):Cou_row(2,nr2),Cou_col(1,nc):Cou_col(2,nc))=Qj2_br;
    CST(Cou_row(1,nr3):Cou_row(2,nr3),Cou_col(1,nc):Cou_col(2,nc))=Qj3_br;     
end
end
end

end

%*******************************
% 2. initial node position matrix
%********************************
function [Node_I,r0]=nodematrix(q,Height,Area)
% 1)geometrical parameters of dome
r0=sqrt(2*Area/(3*sqrt(3)*(9*(2-sqrt(3))*(q-1)^2+3*(q-1)+1)));%initial radius of each unit
t31=3*sqrt(2-sqrt(3))*r0*[sqrt(3)/2;-1/2;0];t32=3*sqrt(2-sqrt(3))*r0*[0;1;0];%location vector

% 2)calculate the center point of each unit Node_C(part0 and 1)
nc=1+q*(q-1);%number of center point(just consider part0 and 1)
Node_C=zeros(3,nc);
Node_C(:,1)=zeros(3,1);
if q>1
    for i=0:q-2
        for k=0:q-1
            n0=1+i*q+k+1;
            Node_C(:,n0)=(i+1)*t31+k*t32;
        end
    end
end

% 3)calculate the nodes of each unit
% the basic unit
NbasicT=zeros(3,3);NbasicB=zeros(3,3);
NbasicT(:,1)=r0*[cos(-pi/12);sin(-pi/12);0];NbasicT(:,2)=r0*[cos(7*pi/12);sin(7*pi/12);0];NbasicT(:,3)=r0*[cos(15*pi/12);sin(15*pi/12);0];
NbasicB(:,1)=r0*[cos(13*pi/12);sin(13*pi/12);0];NbasicB(:,2)=r0*[cos(21*pi/12);sin(21*pi/12);0];NbasicB(:,3)=r0*[cos(5*pi/12);sin(5*pi/12);0];
Nbasic=[NbasicT NbasicB];Nbasic(3,1:3)=Height;
%for part0 and 1
nu=1+3*q*(q-1);%number of unit
Node_I=zeros(3,6*nu);
Node_I(:,1:6*nc)=kron(Node_C,ones(1,6))+kron(ones(1,nc),Nbasic);

% for (2,i,j) and (3,i,j)
c12=cos(2*pi/3);s12=sin(2*pi/3);
C312=[c12 s12 0;-s12 c12 0;0 0 1];
c13=cos(4*pi/3);s13=sin(4*pi/3);
C313=[c13 s13 0;-s13 c13 0;0 0 1];
n11=7;n12=6*(1+(q-1)*q);n21=n12+1;n22=n12+6*(q-1)*q;n31=n22+1;n32=n22+6*(q-1)*q;
Node_I(:,n21:n22)=C312'*Node_I(:,n11:n12);
Node_I(:,n31:n32)=C313'*Node_I(:,n11:n12);

end


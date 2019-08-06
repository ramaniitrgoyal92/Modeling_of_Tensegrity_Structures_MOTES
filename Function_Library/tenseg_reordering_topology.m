% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 

% Re-orders the node indeces if they are not ordered correctely.
% This makes the bars connectivity matrix equal to [-I I 0] in the general
% case
% Also adds virtual nodes for any class-k and pinned node

tenseg_struct_scrambled = tenseg_struct;
ICs_scrambled = ICs;
forces_scrambled = forces;

C_b = field_to_var(tenseg_struct_scrambled,'C_b');
C_s = field_to_var(tenseg_struct_scrambled,'C_s');
N0 = field_to_var(ICs_scrambled,'N0');
Nd0 = field_to_var(ICs_scrambled,'Nd0');
W = field_to_var(forces_scrambled,'W');


if exist('tenseg_struct_scrambled.pinned_nodes','var'),
%     pinned_nodes_scr = field_to_var(tenseg_struct_scrambled,'pinned_nodes');
pinned_nodes_scr = tenseg_struct_scrambled.pinned_nodes_scr;
end

nnodes_or = size(N0,2);
nbars = size(C_b,1);
nstrings = size(C_s,1);

sconn_or = zeros(nstrings,2);

for i = 1 : nstrings

   for j = 1 : nnodes_or

       if C_s(i,j) == -1 

           sconn_or(i,1) = j;

       end

       if C_s(i,j) == 1 

           sconn_or(i,2) = j;

       end

   end

end

class = 0;                        % class of the tensegrity structure
nb_nodes = zeros(nnodes_or,1);    % vector of the number of bars connected to each node
n_virt = 0;                       % total number of virtual nodes
condflag_id = zeros(nnodes_or,1); % vector of id flag of n0 nodes (nodes not connected to bars)
nnodes_0=0;                       % number of nodes not connected to bars

for i = 1 : nnodes_or

    for j = 1 : nbars
    
        nb_nodes(i,1) = nb_nodes(i,1) + abs(C_b(j,i));
        
    end
    
    if class < nb_nodes(i,1)
        
        class = nb_nodes(i,1);
        
    end
    
    if nb_nodes(i,1) > 1
       
        n_virt = n_virt + nb_nodes(i,1) - 1;
        
    end
    
    if nb_nodes(i,1) == 0
       
        condflag_id(i,1) = 1;
        nnodes_0 = nnodes_0 + 1;
        
    end
    
end

nnodes = nnodes_or + n_virt; % total number of nodes (2 * nbars + nnodes_0)
Nmat_NG = zeros(3,nnodes);      % total matrix of nodal coordinates NAGASE

% vector of nodal flags: 
% id = 1: phisical node; 
% id = 2: virtual node; 
% id = 3: n0 node
id_node_flags = zeros(nnodes,1);
id_node_corr = zeros(3,nnodes);  % correspondance between original and new nodes

count = 0;
count_2 = nnodes;
for i = 1 : nnodes_or
    
    if nb_nodes(i,1) == 0
                
        id_node_flags(count_2) = 3;
        
        Nmat_NG(:,count_2) = N0(:,i);
        
        id_node_corr(1,count_2) = i;
        id_node_corr(2,count_2) = count_2;
        
        count_2 = count_2 - 1; 
        
    end
   
    for j = 1 : nb_nodes(i,1)
       
        if j == 1
            
            count = count + 1;
           
            id_node_flags(count) = 1;
            
            Nmat_NG(:,count) = N0(:,i);
            
            id_node_corr(1,count) = i;
            id_node_corr(2,count) = count;
            
        elseif j >= 1
            
            count = count + 1;
            
            id_node_flags(count) = 2;
            
            Nmat_NG(:,count) = N0(:,i);
            
            id_node_corr(1,count) = i;
            id_node_corr(2,count) = count;
                       
        end
        
    end
    
end

id_node_flags = id_node_flags';

CBTmat_or = C_b';
CBmat = zeros(nbars,nnodes);
CBTmat = CBmat';

count = 0;
for i = 1 : nnodes_or
                
    if nb_nodes(i,1) == 1 % one bar connected to ith node

        count = count + 1;

        CBTmat(count,:) = CBTmat_or(i,:);

    end
    
    if nb_nodes(i,1) > 1
                    
        for j = 1 : nbars

            if (abs(CBTmat_or(i,j)) > 0)

                count = count + 1;
                CBTmat(count,j) = CBTmat_or(i,j);

            end

        end
                   
    end
                
end

CBmat_NG = CBTmat';

% reorder connectivity to match Skelton notation

count_1 = 0;
count_2 = nbars;
for i = 1 : nbars
    
    for j = 1 : nnodes

        if CBmat_NG(i,j) == -1

            count_1 = count_1 + 1;

            id_node_corr(3,j) = count_1;

        end
        
        if CBmat_NG(i,j) == +1

            count_2 = count_2 + 1;

            id_node_corr(3,j) = count_2;

        end
    
    end
    
end

for i = 2 * nbars + 1 : nnodes
   
    id_node_corr(3,i) = id_node_corr(2,i);
    
end

Nmat = zeros(3,nnodes);

for i = 1 : nnodes
    
    Nmat(:,id_node_corr(3,i)) = Nmat_NG(:,id_node_corr(2,i));
    
end
ICs.N0 = Nmat;

CBmat = [-eye(nbars) eye(nbars) zeros(nbars,nnodes_0)];
tenseg_struct.C_b = CBmat;

bconn = zeros(nbars,2);

   for i = 1 : nbars
       
       for j = 1 : nnodes
          
           if CBmat(i,j) == 1 
               
               bconn(i,1) = j;
               
           end
           
           if CBmat(i,j) == -1 
               
               bconn(i,2) = j; 
               
           end
           
       end
       
   end
   
sconn = zeros(nstrings,2);

flag_j = 0;
flag_k = 0;
j_new = 0;
k_new = 0;
for i = 1 : nstrings
   
    j_or = sconn_or(i,1);
    k_or = sconn_or(i,2);

    for j = 1 : nnodes
    
        if j_or == id_node_corr(1,j) && flag_j == 0
        
            j_new = id_node_corr(3,j);
            flag_j = 1;
        
        end
             
    end
    
    for k = 1 : nnodes
    
        if (k_or == id_node_corr(1,k)) && (flag_k == 0) && (id_node_corr(3,k) ~= j_new)
        
            k_new = id_node_corr(3,k);
            flag_k = 1;
        
        end
             
    end
    
    sconn(i,1) = j_new;
    sconn(i,2) = k_new;
         
    flag_j = 0;
    flag_k = 0;
    
end

CSTmat = zeros(nnodes,nstrings);
for i = 1 : nstrings
    s1 = sconn(i,1);
    s2 = sconn(i,2);
    CSTmat(s1,i) = 1;
    CSTmat(s2,i) = -1;
end
CSmat = CSTmat';

tenseg_struct.C_s = CSmat;

Crmat = 1/2 * abs(CBmat);
tenseg_struct.C_r = Crmat;

% apply internal constraints (virtual nodes)
if exist('tenseg_struct_scrambled.class_2_nodes','var')
    class_2_nodes = zeros(n_virt,2);
    count=0;
    
    for i = 1 : nnodes
        
        if id_node_flags(1,i) == 2
            
            count = count + 1;
            
            i_node = min(id_node_corr(3,i-1),id_node_corr(3,i));
            j_node = max(id_node_corr(3,i-1),id_node_corr(3,i));
            
            class_2_nodes(count,1) = i_node;
            class_2_nodes(count,2) = j_node;
            
        end
        
    end
    
    tenseg_struct.class_2_nodes = class_2_nodes;
end
    
if exist('tenseg_struct_scrambled.pinned_nodes','var'),
    npin_scr = size(pinned_nodes_scr,1);
    
    pinned_nodes = zeros(1,npin_scr);
    node_flag = 0;
    for i = 1 : npin_scr
        
        for j = 1 : nnodes
            
            if id_node_corr(1,j) == i && node_flag == 0;
                
                node_flag = 1;
                
                pinned_nodes(1,id_node_corr(3,j)) = pinned_nodes_scr(1,i);
                
            end
            
        end
        
        node_flag = 0;
        
    end
end

% compute updated load time hystory

Wmat = zeros(3,nnodes);
node_flag = 0;
for i = 1 : nnodes_or
    
    for j = 1 : nnodes
        
        if id_node_corr(1,j) == i && node_flag == 0;
            
            node_flag = 1;
            
            Wmat(:,id_node_corr(3,j)) = W(:,i);
            
        end
        
    end
    
    node_flag = 0;
    
end

forces.W = Wmat;

% vectors of initial positions and velocities

Ndot_0 = zeros(3,nnodes);

for i = 1 : nnodes
    
    Ndot_0(:,id_node_corr(3,i)) = Nd0(:,id_node_corr(1,i));
    
end

ICs.Nd0 = Ndot_0;

tenseg_struct.m = tenseg_struct_scrambled.m;
tenseg_struct.bar_len_const_hat = tenseg_struct_scrambled.bar_len_const_hat;
if FLAGS.open_loop
    forces.s_0_const = forces_scrambled.s_0_const;
    forces.k = forces_scrambled.k;
    forces.c = forces_scrambled.c;
end





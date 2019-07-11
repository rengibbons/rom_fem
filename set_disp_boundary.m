
function [bc_dofs,disp_dofs,disp_nodes] = set_disp_boundary(sys_vars,mesh_input_file,define_BC,fixed_dof)
%DISP_BOUNDARY_SET Summary of this function goes here
%   Gibbons 12/15/17: This function contains CDISP functionality.

n_dim   = sys_vars.n_dim;
physics = sys_vars.physics;

%% OUTPUT VARIABLEs
bc_dofs      = [];
disp_dofs    = [];
disp_dofs_e  = [0 0];
disp_dofs_c  = [];
disp_nodes   = [];
disp_nodes_e = [];
disp_nodes_c = [];

x_dof_edisp = [];
y_dof_edisp = [];
z_dof_edisp = [];
T_dof_edisp = [];
x_dof_cdisp = [];
y_dof_cdisp = [];
z_dof_cdisp = [];
T_dof_cdisp = [];
x_mag = [];
y_mag = [];
z_mag = [];
T_mag = [];

%%
nodes = sprintf('input/%s.nodes',mesh_input_file);
data  = dlmread(nodes);

%% EDISP implementation
if isfield(define_BC,'EDISP')

    EDISP = define_BC.EDISP;

    nn = 1;
    for tt = 1 : size(EDISP,1)
        nodes_at_edge1(tt,:) = find(data(:,2+EDISP(tt,1))==EDISP(tt,2));
        if strcmp(physics,'elastic') || strcmp(physics,'elastic_neoh')
            x_dof_edisp(nn,:)   = n_dim*(nodes_at_edge1(tt,:)-1) + 1;
            x_dof_edisp(nn+1,:) = EDISP(tt,3)*ones(1,size(x_dof_edisp,2));
            y_dof_edisp(nn,:)   = n_dim*(nodes_at_edge1(tt,:)-1) + 2;
            y_dof_edisp(nn+1,:) = EDISP(tt,4)*ones(1,size(y_dof_edisp,2));
            if n_dim==3
                z_dof_edisp(nn,:)   = n_dim*(nodes_at_edge1(tt,:)-1) + 3;
                z_dof_edisp(nn+1,:) = EDISP(tt,5)*ones(1,size(z_dof_edisp,2));
            end
        elseif strcmp(physics,'thermo_elastic')
            x_dof_edisp(nn,:)   = (n_dim+1)*(nodes_at_edge1(tt,:)-1) + 1;
            x_dof_edisp(nn+1,:) = EDISP(tt,3)*ones(1,size(x_dof_edisp,2));
            y_dof_edisp(nn,:)   = (n_dim+1)*(nodes_at_edge1(tt,:)-1) + 2;
            y_dof_edisp(nn+1,:) = EDISP(tt,4)*ones(1,size(y_dof_edisp,2));
            if n_dim==2
                T_dof_edisp(nn,:)   = (n_dim+1)*(nodes_at_edge1(tt,:)-1) + 3;
                T_dof_edisp(nn+1,:) = EDISP(tt,5)*ones(1,size(T_dof_edisp,2));
            elseif n_dim==3
                z_dof_edisp(nn,:)   = (n_dim+1)*(nodes_at_edge1(tt,:)-1) + 3;
                z_dof_edisp(nn+1,:) = EDISP(tt,5)*ones(1,size(z_dof_edisp,2));
                T_dof_edisp(nn,:)   = (n_dim+1)*(nodes_at_edge1(tt,:)-1) + 4;
                T_dof_edisp(nn+1,:) = EDISP(tt,6)*ones(1,size(T_dof_edisp,2));
            end
        end
        nn = nn + 2;
    end

    if strcmp(physics,'elastic') || strcmp(physics,'elastic_neoh')
        if n_dim==2
            bc_disp_amount = repmat([EDISP(1,3) EDISP(1,4)],size(nodes_at_edge1,2),1);
        end
        if n_dim==3
            bc_disp_amount = repmat([EDISP(1,3) EDISP(1,4) EDISP(1,5)],size(nodes_at_edge1,2),1);
        end
    elseif strcmp(physics,'thermo_elastic')
        if n_dim==2
            bc_disp_amount = repmat([EDISP(1,3) EDISP(1,4) EDISP(1,5)],size(nodes_at_edge1,2),1);
        end
        if n_dim==3
            bc_disp_amount = repmat([EDISP(1,3) EDISP(1,4) EDISP(1,5) EDISP(1,6)],size(nodes_at_edge1,2),1);
        end
    end
    
    bc_dofs_x = reshape(x_dof_edisp,2,[])';
    bc_dofs_y = reshape(y_dof_edisp,2,[])';
    if n_dim==3
        bc_dofs_z = reshape(z_dof_edisp,2,[])';
    end
    if strcmp(physics,'thermo_elastic')
        bc_dofs_T = reshape(T_dof_edisp,2,[])';
    end

    chk1 = isempty(bc_dofs_x);
    chk2 = isempty(bc_dofs_y);
    if n_dim==3
        chk3 = isempty(bc_dofs_z);
    end
    if strcmp(physics,'thermo_elastic')
        chk4 = isempty(bc_dofs_T);
    end
    
    if ~chk1
        disp_dofs_e = union(disp_dofs_e,bc_dofs_x,'rows');
    end
    if ~chk2
        disp_dofs_e = union(disp_dofs_e,bc_dofs_y,'rows');
    end
    if n_dim==3 && ~chk3
        disp_dofs_e = union(disp_dofs_e,bc_dofs_z,'rows');
    end
    if strcmp(physics,'thermo_elastic') && ~chk4
        disp_dofs_e = union(disp_dofs_e,bc_dofs_T,'rows');
    end

    list1 = find(disp_dofs_e(:,2)==0);
    disp_dofs_e(list1,:) = [];

    % OUTPUT
    disp_nodes_e = horzcat(nodes_at_edge1',bc_disp_amount);
%     disp_dofs_e = bc_dofs;

    %% EDISP
%     fixed_dof(ismember(fixed_dof(:,1),bc_dofs(:,1)),:) = [];
%     bc_dofs = union(bc_dofs,fixed_dof,'rows');
end

%% CDSIP implementation
if isfield(define_BC,'CDISP')

    CDISP = define_BC.CDISP;
    disp_nodes_c = zeros(size(CDISP,1),5);
    
    for tt = 1 : size(CDISP,1)
        node_at_specified_coor(tt,1) = find(ismember(data(:,3:3+n_dim-1),CDISP(tt,1:n_dim),'rows'));
        if strcmp(physics,'elastic') || strcmp(physics,'elastic_neoh')
            if CDISP(tt,n_dim+1)
                x_dof_cdisp(tt,:) = n_dim*(node_at_specified_coor(tt,1)-1) + 1;
                x_mag(tt,:)       = CDISP(tt,n_dim+1);
            end
            if CDISP(tt,n_dim+2)
                y_dof_cdisp(tt,:) = n_dim*(node_at_specified_coor(tt,1)-1) + 2;
                y_mag(tt,:)       = CDISP(tt,n_dim+2);
            end
            if n_dim==3 && CDISP(tt,n_dim+3)
                z_dof_cdisp(tt,:) = n_dim*(node_at_specified_coor(tt,1)-1) + 3;
                z_mag(tt,:)       = CDISP(tt,n_dim+3);
            end
        elseif strcmp(physics,'thermo_elastic')
             if CDISP(tt,n_dim+1)
                x_dof_cdisp(tt,:) = (n_dim+1)*(node_at_specified_coor(tt,1)-1) + 1;
                x_mag(tt,:)       = CDISP(tt,n_dim+1);
            end
            if CDISP(tt,n_dim+2)
                y_dof_cdisp(tt,:) = (n_dim+1)*(node_at_specified_coor(tt,1)-1) + 2;
                y_mag(tt,:)       = CDISP(tt,n_dim+2);
            end
            if n_dim==2 && CDISP(tt,n_dim+3)
                T_dof_cdisp(tt,:) = (n_dim+1)*(node_at_specified_coor(tt,1)-1) + 3;
                T_mag(tt,:)       = CDISP(tt,n_dim+3);
            end
            if n_dim==3 && CDISP(tt,n_dim+3)
                z_dof_cdisp(tt,:) = (n_dim+1)*(node_at_specified_coor(tt,1)-1) + 3;
                z_mag(tt,:)       = CDISP(tt,n_dim+3);
            end
            if n_dim==3 && CDISP(tt,n_dim+4)
                T_dof_cdisp(tt,:) = (n_dim+1)*(node_at_specified_coor(tt,1)-1) + 4;
                T_mag(tt,:)       = CDISP(tt,n_dim+4);
            end
        end
        
        disp_nodes_c(tt,1) = node_at_specified_coor(tt);
        if strcmp(physics,'elastic') || strcmp(physics,'elastic_neoh')
            if size(x_mag,1) == tt
                disp_nodes_c(tt,2) = x_mag(tt);
            end
            if size(y_mag,1) == tt
                disp_nodes_c(tt,3) = y_mag(tt);
            end
            if size(z_mag,1) == tt
                disp_nodes_c(tt,4) = z_mag(tt);
            end
        elseif strcmp(physics,'thermo_elastic')
            if size(x_mag,1) == tt
                disp_nodes_c(tt,2) = x_mag(tt);
            end
            if size(y_mag,1) == tt
                disp_nodes_c(tt,3) = y_mag(tt);
            end
            if size(z_mag,1) == tt
                disp_nodes_c(tt,4) = z_mag(tt);
            end
            if size(T_mag,1) == tt
                disp_nodes_c(tt,5) = T_mag(tt);
            end

        end
    end

    if (strcmp(physics,'elastic')  || strcmp(physics,'elastic_neoh')) && n_dim==2
        disp_nodes_c(:,4:5) = [];
    elseif (strcmp(physics,'elastic')  || strcmp(physics,'elastic_neoh')) && n_dim==3
        disp_nodes_c(:,5) = [];
    elseif strcmp(physics,'thermo_elastic') && n_dim==2
        disp_nodes_c(:,4) = [];
    end

%     bc_dofs(find(ismember(bc_dofs,[0. 0.],'rows')),:) = [];

    chk4 = isempty(x_dof_cdisp);
    chk5 = isempty(y_dof_cdisp);
    if n_dim==3
        chk6 = isempty(z_dof_cdisp); 
    end
    chk7 = isempty(T_dof_cdisp);
    
    if n_dim == 2
        fixed_dof_cboun = zeros(size(x_dof_cdisp,1)+size(y_dof_cdisp,1)+size(T_dof_cdisp,1),1);
        fixed_mag_cboun = zeros(size(x_dof_cdisp,1)+size(y_dof_cdisp,1)+size(T_dof_cdisp,1),1);
    elseif n_dim == 3
        fixed_dof_cboun = zeros(size(x_dof_cdisp,1)+size(y_dof_cdisp,1)+size(z_dof_cdisp,1)+size(T_dof_cdisp,1),1);
        fixed_mag_cboun = zeros(size(x_dof_cdisp,1)+size(y_dof_cdisp,1)+size(z_dof_cdisp,1)+size(T_dof_cdisp,1),1);
    end

    start = 1;
    stop  = 1;
    if ~chk4
        stop = size(x_dof_cdisp,1);
        fixed_dof_cboun(start:stop,1) = x_dof_cdisp;
        fixed_mag_cboun(start:stop,1) = x_mag;
    end
    if ~chk5
        start = stop  + 1;
        stop  = start + size(y_dof_cdisp,1) - 1;
        fixed_dof_cboun(start:stop,1) = y_dof_cdisp;
        fixed_mag_cboun(start:stop,1) = y_mag;
    end
    if n_dim==3 && ~chk6
        start = stop  + 1;
        stop  = start + size(y_dof_cdisp,1) - 1;
        fixed_dof_cboun(start:stop,1) = z_dof_cdisp;
        fixed_mag_cboun(start:stop,1) = z_mag;
    end
    if strcmp(physics,'thermo_elastic') && ~chk7
        start = stop  + 1;
        stop  = start + size(y_dof_cdisp,1) - 1;
        fixed_dof_cboun(start:stop,1) = T_dof_cdisp;
        fixed_mag_cboun(start:stop,1) = T_mag;
    end
    
    while fixed_dof_cboun(1,1)==0
        fixed_dof_cboun(1,:) = [];
        fixed_mag_cboun(1,:) = [];
    end
    
    disp_dofs_c = sortrows([fixed_dof_cboun,fixed_mag_cboun]);

end

%% OUTPUT

disp_nodes = vertcat(disp_nodes_e,disp_nodes_c);
disp_dofs  = vertcat(disp_dofs_e,disp_dofs_c);

% CDISP
bc_dofs = union(fixed_dof,disp_dofs,'rows');

while bc_dofs(1,1)==0
    bc_dofs(1,:) = [];
end

end













% %% OUTPUT VARIABLEs
% bc_dofs    = [0 0];
% disp_dofs  = [];
% disp_nodes = [];
% 
% x_dof_cboun = [];
% y_dof_cboun = [];
% z_dof_cboun = [];
% T_dof_cboun = [];
% x_mag = [];
% y_mag = [];
% z_mag = [];
% T_mag = [];
% %%
% fixed_dofs      = zeros(length(fixed_dof),2);
% fixed_dofs(:,1) = fixed_dof;
% 
% nodes = sprintf('input/%s.nodes',mesh_input_file);
% data  = dlmread(nodes);
% 
% %% EDISP implementation
% if isfield(define_BC,'EDISP') || isfield(define_BC,'CDISP')
%     
%     if isfield(define_BC,'EDISP')
%         
%         EDISP = define_BC.EDISP; %[2 0.5 0 1];
%             
%         nn = 1;
%         for tt = 1 : size(EDISP,1)
%             nodes_at_edge1(tt,:) = find(data(:,2+EDISP(tt,1))==EDISP(tt,2));
%             x_dof(nn,:)   = n_dim*(nodes_at_edge1(tt,:)-1) + 1;
%             x_dof(nn+1,:) = zeros(1,size(x_dof,2))        + EDISP(tt,3);
%             y_dof(nn,:)   = n_dim*(nodes_at_edge1(tt,:)-1) +2;
%             y_dof(nn+1,:) = zeros(1,size(y_dof,2))        + EDISP(tt,4);
% 
%             if n_dim==3
%                 z_dof(nn,:)   = n_dim*(nodes_at_edge1(tt,:)-1) + 3;
%                 z_dof(nn+1,:) = zeros(1,size(z_dof,2))        + EDISP(tt,5);
%             end
%             nn = nn + 2;
%         end
%         
%         if n_dim==2 
%             bc_disp_amount = repmat([EDISP(1,3) EDISP(1,4)],size(nodes_at_edge1,2),1);
%         end
%         if n_dim==3
%             bc_disp_amount = repmat([EDISP(1,3) EDISP(1,4) EDISP(1,5)],size(nodes_at_edge1,2),1);
%         end
% 
%         bc_dofs_x = reshape(x_dof,2,[])';
%         bc_dofs_y = reshape(y_dof,2,[])';
%         if n_dim==3
%             bc_dofs_z = reshape(z_dof,2,[])';
%         end
%         
%         chk1 = isempty(bc_dofs_x);
%         chk2 = isempty(bc_dofs_y);
%         if n_dim==3
%             chk3 = isempty(bc_dofs_z);
%         end
% 
%         if ~chk1
%             bc_dofs = union(bc_dofs,bc_dofs_x,'rows');
%         end
%         if ~chk2
%             bc_dofs = union(bc_dofs,bc_dofs_y,'rows');
%         end
%         if n_dim==3 && ~chk3
%             bc_dofs = union(bc_dofs,bc_dofs_z,'rows');
%         end
%         
%         list1 = find(bc_dofs(:,2)==0);
%         bc_dofs(list1,:) = [];
%         
%         % OUTPUT
%         disp_nodes = horzcat(nodes_at_edge1',bc_disp_amount);
%         
%         % OUTPUT
%         disp_dofs = bc_dofs(:,1);
%         
%         %% EDISP
%         fixed_dofs(ismember(fixed_dofs(:,1),bc_dofs(:,1)),:) = [];
%         
%         bc_dofs = union(bc_dofs,fixed_dofs,'rows');
%         
%     end
%     
%     if (isfield(define_BC,'CDISP'))
%         
%         CDISP = define_BC.CDISP; %=[2.0 -0.5 1 1];'
%         disp_nodes_c = zeros(size(CDISP,1),4);
%         
%         for tt = 1 : size(CDISP,1)
%             node_at_specified_coor(tt,1) = find(ismember(data(:,3:3+n_dim-1),CDISP(tt,1:n_dim),'rows'));
%             if CDISP(tt,n_dim+1)
%                 x_dof_cboun(tt,:) = n_dim*(node_at_specified_coor(tt,1)-1)+1;
%                 x_mag(tt,:)       = CDISP(tt,n_dim+1);
%             end
%             if CDISP(tt,n_dim+2)
%                 y_dof_cboun(tt,:) = n_dim*(node_at_specified_coor(tt,1)-1)+2;
%                 y_mag(tt,:)       = CDISP(tt,n_dim+2);
%             end
%             if n_dim==3 && CDISP(tt,n_dim+3)
%                 z_dof_cboun(tt,:) = n_dim*(node_at_specified_coor(tt,1)-1)+3;
%                 z_mag(tt,:)       = CDISP(tt,n_dim+3);
%             end
%             
%             disp_nodes_c(tt,1) = node_at_specified_coor(tt);
%             if size(x_mag,1) == tt
%                 disp_nodes_c(tt,2) = x_mag(tt);
%             end
% 
%             if size(y_mag,1) == tt
%                 disp_nodes_c(tt,3) = y_mag(tt);
%             end
% 
%             if size(z_mag,1) == tt
%                 disp_nodes_c(tt,4) = z_mag(tt);
%             end
%             
%             if n_dim == 2
%                 disp_nodes_c = disp_nodes_c(:,1:3);
%             end
%         end
%         
%         bc_dofs(find(ismember(bc_dofs,[0. 0.],'rows')),:) = [];
% 
%         chk4 = isempty(x_dof_cboun);
%         chk5 = isempty(y_dof_cboun);
%         if n_dim==3
%             chk6=isempty(z_dof_cboun); 
%         end
% 
%         if n_dim == 2
%             fixed_dofs_cboun = zeros(size(x_dof_cboun,1)+size(y_dof_cboun,1),1);
%             fixed_mag_cboun  = zeros(size(x_dof_cboun,1)+size(y_dof_cboun,1),1);
%         elseif n_dim == 3
%             fixed_dofs_cboun = zeros(size(x_dof_cboun,1)+size(y_dof_cboun,1)+size(z_dof_cboun,1),1);
%             fixed_mag_cboun  = zeros(size(x_dof_cboun,1)+size(y_dof_cboun,1)+size(z_dof_cboun,1),1);
%         end
%         
%         start = 1;
%         stop  = 1;
%         if ~chk4
%             stop = size(x_dof_cboun,1);
%             fixed_dofs_cboun(start:stop,1) = x_dof_cboun;
%             fixed_mag_cboun(start:stop,1)  = x_mag;
%         end
%         
%         if ~chk5
%             start = stop  + 1;
%             stop  = start + size(y_dof_cboun,1) - 1;
%             fixed_dofs_cboun(start:stop,1) = y_dof_cboun;
%             fixed_mag_cboun(start:stop,1)  = y_mag;
%         end
%         
%         if n_dim==3 && ~chk6
%             start = stop  + 1;
%             stop  = start + size(y_dof_cboun,1) - 1;
%             fixed_dofs_cboun(start:stop,1) = z_dof_cboun;
%             fixed_mag_cboun(start:stop,1)  = y_mag;
%         end
%         
%         if n_dim==2
%             bc_disp_amount = repmat([CDISP(1,3) CDISP(1,4)],length(fixed_dofs_cboun),1); 
%         end
%         if n_dim==3
%             bc_disp_amount = repmat([CDISP(1,3) CDISP(1,4) CDISP(1,5)],length(fixed_dofs_cboun),1);
%         end
%         
%         % OUTPUT
%         disp_nodes = vertcat(disp_nodes,disp_nodes_c);
% 
%         % OUTPUT
%         disp_dofs = [fixed_dofs_cboun,fixed_mag_cboun];
%         
%         %% CDISP        
%         bc_dofs = union(bc_dofs,fixed_dofs,'rows');
%         bc_dofs = union(bc_dofs,disp_dofs,'rows');
%     end
%     
% else
%     bc_dofs = fixed_dofs;
% end
% 
% bc_dofs = bc_dofs(2:end,:);
% 
% end


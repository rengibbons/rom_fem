function [fixed_dofs,fixed_nodes] = set_dirichlet_boundary(sys_vars,mesh_input_file,define_BC,dof)

n_dim   = sys_vars.n_dim;
physics = sys_vars.physics;
if strcmp(physics,'thermo_elastic')
    T0 = sys_vars.mat_para.T0;
end

%% OUTPUT VALUES
fixed_dofs         = [];
fixed_nodes        = [];
fixed_dofs_eboun   = [];
fixed_dofs_eboun_T = [];
fixed_dofs_cboun   = [];
fixed_dofs_cboun_T = [];
fixed_nodes_cboun  = [];
x_dof = [];
y_dof = [];
z_dof = [];
T_dof = [];
x_dof_cboun = [];
y_dof_cboun = [];
z_dof_cboun = [];
T_dof_cboun = [];

%%
nodes = sprintf('input/%s.nodes',mesh_input_file);
data  = dlmread(nodes);

%% EBOUN implementation
if (isfield(define_BC,'EBOUN'))       
    EBOUN = define_BC.EBOUN;

    for tt = 1 : size(EBOUN,1)
        nodes_at_edge1(tt,:) = find(data(:,2+EBOUN(tt,1))==EBOUN(tt,2));
        if strcmp(physics,'elastic') || strcmp(physics,'elastic_neoh')
            if EBOUN(tt,3)
                x_dof(tt,:) = n_dim*(nodes_at_edge1(tt,:)-1) + 1;
            end
            if EBOUN(tt,4)
                y_dof(tt,:) = n_dim*(nodes_at_edge1(tt,:)-1) + 2;
            end
            if n_dim==3 && EBOUN(tt,5)
                z_dof(tt,:) = n_dim*(nodes_at_edge1(tt,:)-1) + 3;
            end
        elseif strcmp(physics,'thermo_elastic')
            if EBOUN(tt,3)
                x_dof(tt,:) = (n_dim+1)*(nodes_at_edge1(tt,:)-1) + 1;
            end
            if EBOUN(tt,4)
                y_dof(tt,:) = (n_dim+1)*(nodes_at_edge1(tt,:)-1) + 2;
            end
            if n_dim==2 && EBOUN(tt,5)
                T_dof(tt,:) = (n_dim+1)*(nodes_at_edge1(tt,:)-1) + 3;
            end
            if n_dim==3 && EBOUN(tt,5)
                z_dof(tt,:) = (n_dim+1)*(nodes_at_edge1(tt,:)-1) + 3;
            end
            if n_dim==3 && EBOUN(tt,6)
                T_dof(tt,:) = (n_dim+1)*(nodes_at_edge1(tt,:)-1) + 4;
            end
        end
    end
    
    chk1 = isempty(x_dof);
    chk2 = isempty(y_dof);
    if n_dim==3
        chk3 = isempty(z_dof);
    end
    chk4 = isempty(T_dof);
    
    if ~chk1
        fixed_dofs_eboun = union(fixed_dofs_eboun,x_dof');
    end
    if ~chk2
        fixed_dofs_eboun = union(fixed_dofs_eboun,y_dof');
    end
    if n_dim==3 && ~chk3
        fixed_dofs_eboun = union(fixed_dofs_eboun,z_dof');
    end
    
    fixed_dofs_eboun = [fixed_dofs_eboun,zeros(size(fixed_dofs_eboun,1),1)];
    
    if ~chk4
        fixed_dofs_eboun_T = [T_dof',T0*ones(size(T_dof',1),1)];
%         fixed_dofs = union(fixed_dofs,T_dof');
    end
    
    fixed_dofs_eboun = sortrows([fixed_dofs_eboun;fixed_dofs_eboun_T]);
    fixed_nodes = unique(sort(reshape(nodes_at_edge1,[],1)));
end

%% CBOUN implementation
if (isfield(define_BC,'CBOUN'))
    CBOUN = define_BC.CBOUN;
    
    for tt = 1 : size(CBOUN,1)
        node_at_specified_coor(tt,1) = find(ismember(data(:,3:3+n_dim-1),CBOUN(tt,1:n_dim),'rows'));
        if strcmp(physics,'elastic') || strcmp(physics,'elastic_neoh')
            if CBOUN(tt,n_dim+1)
                x_dof_cboun(tt,:) = n_dim*(node_at_specified_coor(tt,1)-1) + 1;
            end
            if CBOUN(tt,n_dim+2)
                y_dof_cboun(tt,:) = n_dim*(node_at_specified_coor(tt,1)-1) + 2;
            end
            if n_dim==3 && CBOUN(tt,n_dim+3)
                z_dof_cboun(tt,:) = n_dim*(node_at_specified_coor(tt,1)-1) + 3;
            end
        elseif strcmp(physics,'thermo_elastic')
            if CBOUN(tt,n_dim+1)
                x_dof_cboun(tt,:) = (n_dim+1)*(node_at_specified_coor(tt,1)-1) + 1;
            end
            if CBOUN(tt,n_dim+2)
                y_dof_cboun(tt,:) = (n_dim+1)*(node_at_specified_coor(tt,1)-1) + 2;
            end
            if n_dim==2 && CBOUN(tt,n_dim+3)
                T_dof_cboun(tt,:) = (n_dim+1)*(node_at_specified_coor(tt,1)-1) + 3;
            end
            if n_dim==3 && CBOUN(tt,n_dim+3)
                z_dof_cboun(tt,:) = (n_dim+1)*(node_at_specified_coor(tt,1)-1) + 3;
            end
            if n_dim==3 && CBOUN(tt,n_dim+4)
                T_dof_cboun(tt,:) = (n_dim+1)*(node_at_specified_coor(tt,1)-1) + 4;
            end
        end
    end
    
    chk4 = isempty(x_dof_cboun);
    chk5 = isempty(y_dof_cboun);
    if n_dim==3
        chk6 = isempty(z_dof_cboun);
    end
    chk7 = isempty(T_dof_cboun);
    
    fixed_dofs_cboun = [];

    if ~chk4
%         fixed_dofs_cboun = union(fixed_dofs_cboun,[unique(x_dof_cboun),zeros(size(unique(x_dof_cboun),1),1)]);
        fixed_dofs_cboun = union(fixed_dofs_cboun,x_dof_cboun');
    end
    
    if ~chk5
        fixed_dofs_cboun = union(fixed_dofs_cboun,y_dof_cboun');
%         fixed_dofs_cboun = union(fixed_dofs_cboun,[unique(y_dof_cboun),zeros(size(unique(y_dof_cboun),1),1)]);
    end

    if n_dim==3 && ~chk6
        fixed_dofs_cboun = union(fixed_dofs_cboun,z_dof_cboun');
%         fixed_dofs_cboun = union(fixed_dofs_cboun,[unique(z_dof_cboun),zeros(size(unique(z_dof_cboun),1),1)]);
    end
    
    fx_sz = size(fixed_dofs_cboun);
    if fx_sz(1)==1 && fx_sz(2)>1
        fixed_dofs_cboun = fixed_dofs_cboun';
    end
    fixed_dofs_cboun = [fixed_dofs_cboun,zeros(size(fixed_dofs_cboun,1),1)];

    if ~chk7
        fixed_dofs_cboun_T = [T_dof_cboun,T0*ones(size(T_dof_cboun,1),1)];
%         fixed_dofs_cboun = union(fixed_dofs_cboun,T_dof_cboun');
%         fixed_dofs_cboun = union(fixed_dofs_cboun,[unique(T_dof_cboun),T0*ones(size(unique(x_dof_cboun),1),1)]);
    end

    fixed_dofs_cboun  = sortrows([fixed_dofs_cboun;fixed_dofs_cboun_T]); 
    fixed_nodes_cboun = node_at_specified_coor;

end

% fixed_dofs        = union(fixed_dofs,fixed_dofs_cboun);
fixed_dofs =  sortrows([fixed_dofs_eboun;fixed_dofs_cboun]);

while fixed_dofs(1,1)==0
    fixed_dofs(1,:) = [];
end
fixed_dofs = unique(fixed_dofs,'rows');
fixed_nodes       = union(fixed_nodes,fixed_nodes_cboun);
end

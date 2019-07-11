function [fixed_dofs,fixed_nodes] = dirichlet_boundary_set(mesh_input_file,geometry,dof)

global ndim;

bounds=sprintf('input/%s.bouns',mesh_input_file);
nodes=sprintf('input/%s.nodes',mesh_input_file);

% file1 = fopen(bounds,'w');

data = dlmread(nodes);

    %% CUBE needs modification
if (strcmp(geometry,'cube'))
    
    %% SET BOTTOM
    edge_data=data(:,5);
    % be careful! always different!
    lowest_value_edge_1= min(edge_data)+0.01;
    nodes_at_edge1= find(edge_data<lowest_value_edge_1);
    
    is_x_fixed=1; x_dof=[0 0];
    is_y_fixed=1; y_dof=[0 0];
    is_z_fixed=1; z_dof=[0 0];
    
    if (is_x_fixed); x_dof=ndim*(nodes_at_edge1-1)+1;  end;
    if (is_y_fixed); y_dof=ndim*(nodes_at_edge1-1)+2;  end;
    if (is_z_fixed); z_dof=ndim*(nodes_at_edge1-1)+3;  end;
    
    fixed_dofs_bottom= union(x_dof,y_dof);
    fixed_dofs_bottom= union(fixed_dofs_bottom,z_dof);
    fixed_nodes_bottom=nodes_at_edge1;
    
    %% SET TOP
    edge_data=data(:,5);
    % be careful! always different!
    lowest_value_edge_1= max(edge_data)-0.01;
    nodes_at_edge1= find(edge_data>lowest_value_edge_1);
    
    is_x_fixed=0; x_dof=[0];
    is_y_fixed=0; y_dof=[0];
    is_z_fixed=0; z_dof=[0];
    
    if (is_x_fixed); x_dof=ndim*(nodes_at_edge1-1)+1;  end;
    if (is_y_fixed); y_dof=ndim*(nodes_at_edge1-1)+2;  end;
    if (is_z_fixed); z_dof=ndim*(nodes_at_edge1-1)+3;  end;
    
    fixed_dofs_top= union(x_dof,y_dof);
    fixed_dofs_top= union(fixed_dofs_top,z_dof);
    indices = find(fixed_dofs_top==0);
    fixed_dofs_top(indices) = [];
    fixed_nodes_top=nodes_at_edge1;
    %%
    
    if (isempty(fixed_dofs_top))
        fixed_dofs=fixed_dofs_bottom;
        fixed_nodes=fixed_nodes_bottom;
    else
        fixed_dofs=union(fixed_dofs_top,fixed_dofs_bottom);
        fixed_nodes=union(fixed_nodes_top,fixed_nodes_bottom);
    end
    
    
    %% CIRCULAR
elseif (strcmp(geometry,'circular'))
    
    a= data(:,3).^2+data(:,4).^2 - 60^2;
    bottom= min(data(:,5))+1;
    loca1= find(abs(a)<10^-1);
    loca2= find( data(:,5)<bottom );
    
    bounds = data(loca1,:);
    bottom_bounds = data(loca2,:);
    
    for i=1:length(bounds)
        fprintf(file1,'%d  0  1 1 0\n',bounds(i,1));
    end
    for i=1:length(bottom_bounds)
        fprintf(file1,'%d  0  1 1 1\n',bottom_bounds(i,1));
    end
    fprintf(file1,'\n');
    
    %% QUAD CASE
elseif (strcmp(geometry,'quad'))
    EBOUN=[2 -0.5 1 1];
    
    for tt=1:size(EBOUN,1)
        nodes_at_edge1(tt,:)=find(data(:,2+EBOUN(tt,1))==EBOUN(tt,2));
        if (EBOUN(tt,3)); x_dof(tt,:)=ndim*(nodes_at_edge1(tt,:)-1)+1;  end;
        if (EBOUN(tt,4)); y_dof(tt,:)=ndim*(nodes_at_edge1(tt,:)-1)+2;  end;
    end
    
    chk1=isempty(x_dof);
    chk2=isempty(y_dof);
    
    if (~chk1 && ~chk2)
        fixed_dofs=sort(union(x_dof,y_dof))';
    else
        if (chk1) fixed_dofs=sort(reshape(x_dof,[],1))'; end;
        if (chk2) fixed_dofs=sort(reshape(y_dof,[],1))'; end;
    end
    fixed_nodes=sort(reshape(nodes_at_edge1,[],1));
    %% ERROR
else
    fprintf(file1,'WORNG BOUNDARY CONDITION OPTIONS');
end
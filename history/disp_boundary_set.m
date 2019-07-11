function [bc_dofs,disp_dofs,disp_nodes] = disp_boundary_set(mesh_input_file,geometry,dof_list,fixed_dof)
%DISP_BOUNDARY_SET Summary of this function goes here
%   Detailed explanation goes here
global ndim;

bounds=sprintf('input/%s.bouns',mesh_input_file);
nodes=sprintf('input/%s.nodes',mesh_input_file);

% file1 = fopen(bounds,'w');

data = dlmread(nodes);

if (strcmp(geometry,'cube'))
    
    edge_data=data(:,5);
    % be careful! always different!
    lowest_value_edge_1= max(edge_data)-0.01;
    nodes_at_edge1= find(edge_data>lowest_value_edge_1);
    num_nodes_at_edge1=length(nodes_at_edge1);
    disp_vector=zeros(num_nodes_at_edge1,1);
    %% SET THE DISPLACEMENT VALUES AT DOFS (X,Y,Z)
    is_x_disp=0;       x_dof=[0 0]; disp_x=0;
    is_y_disp=0;       y_dof=[0 0]; disp_y=0;
    is_z_disp=1;       z_dof=[0 0]; disp_z=0.5;
    %%
    if (is_x_disp); x_dof=[ndim*(nodes_at_edge1-1)+1 , disp_vector+disp_x ];  end;
    if (is_y_disp); y_dof=[ndim*(nodes_at_edge1-1)+2 , disp_vector+disp_y ];  end;
    if (is_z_disp); z_dof=[ndim*(nodes_at_edge1-1)+3 , disp_vector+disp_z ];  end;
    
    disp_dofs= union(x_dof,y_dof,'rows');
    disp_dofs= union(disp_dofs,z_dof,'rows');
    indices = find(disp_dofs(:,1)==0);
    disp_dofs(indices,:) = [];
    
    disp_nodes=nodes_at_edge1;
    
    fixed_dofs=zeros(length(fixed_dof),2);
    fixed_dofs(:,1)=fixed_dof;
    
    if (~isempty(disp_dofs))
        bc_dofs=union(disp_dofs,fixed_dofs,'rows');
    else
        disp_nodes=[];
        bc_dofs=fixed_dofs;
    end
    %% QUAD
elseif (strcmp(geometry,'quad'))
    %            coor location amount dof1 amount dof2
    EDISP=[2 0.5 0 1];
    
    nn=1;
    for tt=1:size(EDISP,1)
        nodes_at_edge1(tt,:)=find(data(:,2+EDISP(tt,1))==EDISP(tt,2));
        x_dof(nn,:)=ndim*(nodes_at_edge1(tt,:)-1)+1;
        x_dof(nn+1,:)=zeros(1,size(x_dof,2))+EDISP(tt,3);
        y_dof(nn,:)=ndim*(nodes_at_edge1(tt,:)-1)+2;
        y_dof(nn+1,:)=zeros(1,size(y_dof,2))+EDISP(tt,4);
        nn=nn+2;
    end
    
    bc_dofs_x=reshape(x_dof,2,[])';
    bc_dofs_y=reshape(y_dof,2,[])';
    
    chk1=isempty(bc_dofs_x);
    chk2=isempty(bc_dofs_y);
    
    if (~chk1 && ~chk2)
        bc_dofs=union(bc_dofs_x,bc_dofs_y,'rows');
    else
        if (chk1) bc_dofs=bc_dofs_x; end;
        if (chk2) bc_dofs=bc_dofs_y; end;
    end
    
    list1=find(bc_dofs(:,2)==0);
    bc_dofs(list1,:)=[];
    disp_dofs=bc_dofs;
    disp_nodes=sort(reshape(nodes_at_edge1,[],1));
    
    %% ERROR
end

    
    fixed_dofs=zeros(length(fixed_dof),2);
    fixed_dofs(:,1)=fixed_dof;


    if (~isempty(disp_dofs))
        bc_dofs=union(disp_dofs,fixed_dofs,'rows');
    else
        disp_nodes=[];
        bc_dofs=fixed_dofs;
    end

end


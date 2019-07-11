function [sys_vars,coor_list,connect_list,element_dof,dof] = mesh(sys_vars,mesh_input_file)

% (MESH + BOUNDARY C. + LOADING) DESCRIPTIONS + GRAPHICAL OUTPUT
file01  = sys_vars.file01;
physics = sys_vars.physics;
if strcmp(physics,'thermo_elastic')
    T0 = sys_vars.mat_para.T0;
end

% 1) GET THE MESH AS FEAP GET FROM CUBIT
nodes = sprintf('input/%s.nodes',mesh_input_file);
elems = sprintf('input/%s.elems',mesh_input_file);

% FOR 3D MESH, READ NODE COORDINATES
coor_list    = dlmread(nodes);
connect_list = dlmread(elems);

total_element_no = size(connect_list,1);
total_node_no    = size(coor_list,1);
node_p_elm       = size(connect_list,2)-3;

if strcmp(physics,'elastic') || strcmp(physics,'elastic_neoh')
    dof         = zeros(total_node_no*sys_vars.n_dim,1);
    element_dof = zeros(total_element_no,node_p_elm*sys_vars.n_dim);
    
    if (sys_vars.n_dim==3)
        x=sys_vars.n_dim*(connect_list(:,4:end)-1)+1;
        y=sys_vars.n_dim*(connect_list(:,4:end)-1)+2;
        z=sys_vars.n_dim*(connect_list(:,4:end)-1)+3;

        element_dof(:,1:3:end)=x(:,:);
        element_dof(:,2:3:end)=y(:,:);
        element_dof(:,3:3:end)=z(:,:);

        dof(1:sys_vars.n_dim:end,1)=coor_list(:,3);
        dof(2:sys_vars.n_dim:end,1)=coor_list(:,4);
        dof(3:sys_vars.n_dim:end,1)=coor_list(:,5);
        
    elseif (sys_vars.n_dim==2)
        x=sys_vars.n_dim*(connect_list(:,4:end)-1)+1;
        y=sys_vars.n_dim*(connect_list(:,4:end)-1)+2;

        element_dof(:,1:2:end)=x(:,:);
        element_dof(:,2:2:end)=y(:,:);

        dof(1:sys_vars.n_dim:end,1)=coor_list(:,3);
        dof(2:sys_vars.n_dim:end,1)=coor_list(:,4);
        
    end  
     
elseif strcmp(physics,'thermo_elastic')
    dof         = zeros(total_node_no*(sys_vars.n_dim+1),1);
    element_dof = zeros(total_element_no,node_p_elm*(sys_vars.n_dim+1));
    coor_list  = horzcat(coor_list,T0*ones(size(coor_list,1),1));
    
    if (sys_vars.n_dim==3)
        % Untested, but should work. RG 12/30/17
        x = (sys_vars.n_dim+1)*(connect_list(:,4:end)-1)+1;
        y = (sys_vars.n_dim+1)*(connect_list(:,4:end)-1)+2;
        z = (sys_vars.n_dim+1)*(connect_list(:,4:end)-1)+3;
        T = (sys_vars.n_dim+1)*(connect_list(:,4:end)-1)+4;

        element_dof(:,1:4:end) = x(:,:);
        element_dof(:,2:4:end) = y(:,:);
        element_dof(:,3:4:end) = z(:,:);
        element_dof(:,4:4:end) = T(:,:);

        dof(1:sys_vars.n_dim+1:end,1) = coor_list(:,3);
        dof(2:sys_vars.n_dim+1:end,1) = coor_list(:,4);
        dof(3:sys_vars.n_dim+1:end,1) = coor_list(:,5);
        dof(4:sys_vars.n_dim+1:end,1) = T0;
        
    elseif (sys_vars.n_dim==2)
        x = (sys_vars.n_dim+1)*(connect_list(:,4:end)-1)+1;
        y = (sys_vars.n_dim+1)*(connect_list(:,4:end)-1)+2;
        T = (sys_vars.n_dim+1)*(connect_list(:,4:end)-1)+3;
        
        element_dof(:,1:3:end) = x(:,:);
        element_dof(:,2:3:end) = y(:,:);
        element_dof(:,3:3:end) = T(:,:);

        dof(1:sys_vars.n_dim+1:end,1) = coor_list(:,3);
        dof(2:sys_vars.n_dim+1:end,1) = coor_list(:,4);
        dof(3:sys_vars.n_dim+1:end,1) = T0;
        
    end  
end

% COMPATIBILITY FOR FEAP
% coor_list = [NODE,0, X   Y   Z]
% connect_list = [ELEMENT,0,MAT_NO,CONNECTIVITY OF NODES]
% fprintf(file01,'\nElement type: %10s \t Node Number per element: %2d \nTotal node number: %6d \t Total element number: %6d\n\n ',sys_vars.element_type,node_p_elm,total_node_no,total_element_no);

sys_vars.total_element_no = total_element_no;
sys_vars.total_node_no    = total_node_no;
sys_vars.node_p_elm       = node_p_elm;

end

%% This function works for writing the mesh with no thermoelastic functionality.
% function [sys_vars,coor_list,connect_list,element_dof,dof] = mesh(sys_vars,mesh_input_file)
% 
% % (MESH + BOUNDARY C. + LOADING) DESCRIPTIONS + GRAPHICAL OUTPUT
% file01 = sys_vars.file01;
% 
% 
% % 1) GET THE MESH AS FEAP GET FROM CUBIT
% nodes = sprintf('input/%s.nodes',mesh_input_file);
% elems = sprintf('input/%s.elems',mesh_input_file);
% 
% % FOR 3D MESH, READ NODE COORDINATES
% coor_list    = dlmread(nodes);
% connect_list = dlmread(elems);
% 
% total_element_no = size(connect_list,1);
% total_node_no    = size(coor_list,1);
% node_p_elm       = size(connect_list,2)-3;
% 
% dof         = zeros(total_node_no*sys_vars.n_dim,1);
% element_dof = zeros(total_element_no,node_p_elm*sys_vars.n_dim);
% 
% if (sys_vars.n_dim==3)
%     
%     x=sys_vars.n_dim*(connect_list(:,4:end)-1)+1;
%     y=sys_vars.n_dim*(connect_list(:,4:end)-1)+2;
%     z=sys_vars.n_dim*(connect_list(:,4:end)-1)+3;
%     
%     element_dof(:,1:3:end)=x(:,:);
%     element_dof(:,2:3:end)=y(:,:);
%     element_dof(:,3:3:end)=z(:,:);
%     
%     dof(1:sys_vars.n_dim:end,1)=coor_list(:,3);
%     dof(2:sys_vars.n_dim:end,1)=coor_list(:,4);
%     dof(3:sys_vars.n_dim:end,1)=coor_list(:,5);
%     
% elseif (sys_vars.n_dim==2)
%         
%     x=sys_vars.n_dim*(connect_list(:,4:end)-1)+1;
%     y=sys_vars.n_dim*(connect_list(:,4:end)-1)+2;
%     
%     element_dof(:,1:2:end)=x(:,:);
%     element_dof(:,2:2:end)=y(:,:);
%     
%     dof(1:sys_vars.n_dim:end,1)=coor_list(:,3);
%     dof(2:sys_vars.n_dim:end,1)=coor_list(:,4);
%     
% end
% 
% % coor_list_check(:,1)=dof(1:sys_vars.n_dim:end,1)
% % coor_list_check(:,2)=dof(2:sys_vars.n_dim:end,1)
% % coor_list_check(:,3)=dof(3:sys_vars.n_dim:end,1)
% 
% % COMPATIBILITY FOR FEAP
% % coor_list = [NODE,0, X   Y   Z]
% % connect_list = [ELEMENT,0,MAT_NO,CONNECTIVITY OF NODES]
% fprintf(file01,'\nElement type: %10s \t Node Number per element: %2d \nTotal node number: %6d \t Total element number: %6d\n\n ',sys_vars.element_type,node_p_elm,total_node_no,total_element_no);
% 
% sys_vars.total_element_no = total_element_no;
% sys_vars.total_node_no    = total_node_no;
% sys_vars.node_p_elm       = node_p_elm;
%  
% end
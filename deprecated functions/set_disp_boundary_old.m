function [bc_dofs,disp_dofs,disp_nodes] = set_disp_boundary_old(ndim,mesh_input_file,define_BC,dof_list,fixed_dof)
%DISP_BOUNDARY_SET Summary of this function goes here
%   Detailed explanation goes here
%   Gibbons 12/15/17: This function does not contain CDISP functionality.
%                     Use set_disp_boundary_2() instead.

%% OUTPUT VARIABLEs
bc_dofs    = [0 0];
disp_dofs  = [];
disp_nodes = [];

%%
fixed_dofs      = zeros(length(fixed_dof),2);
fixed_dofs(:,1) = fixed_dof;

nodes = sprintf('input/%s.nodes',mesh_input_file);
data  = dlmread(nodes);

%% EDISP implementation
if (isfield(define_BC,'EDISP'))
    
    EDISP = define_BC.EDISP; %[2 0.5 0 1];
    
    nn = 1;
    for tt = 1 : size(EDISP,1)
        nodes_at_edge1(tt,:) = find(data(:,2+EDISP(tt,1))==EDISP(tt,2));
        x_dof(nn,:)          = ndim*(nodes_at_edge1(tt,:)-1)+1;
        x_dof(nn+1,:)        = zeros(1,size(x_dof,2))+EDISP(tt,3);
        y_dof(nn,:)          = ndim*(nodes_at_edge1(tt,:)-1)+2;
        y_dof(nn+1,:)        = zeros(1,size(y_dof,2))+EDISP(tt,4);
        %
        if(ndim==3)
            z_dof(nn,:)   = ndim*(nodes_at_edge1(tt,:)-1)+3;
            z_dof(nn+1,:) = zeros(1,size(z_dof,2))+EDISP(tt,5);
        end
        %
        nn = nn + 2;
    end
    
    if(ndim==2)
        bc_disp_amount = repmat([EDISP(1,3) EDISP(1,4)],size(nodes_at_edge1,2),1); 
    end
    if(ndim==3)
        bc_disp_amount = repmat([EDISP(1,3) EDISP(1,4) EDISP(1,5)],size(nodes_at_edge1,2),1);
    end
   
    bc_dofs_x = reshape(x_dof,2,[])';
    bc_dofs_y = reshape(y_dof,2,[])';
    if (ndim==3)
        bc_dofs_z=reshape(z_dof,2,[])';
    end
    
    chk1 = isempty(bc_dofs_x);
    chk2 = isempty(bc_dofs_y);
    if (ndim==3)
        chk3 = isempty(bc_dofs_z);
    end;
    if (~chk1)
        bc_dofs = union(bc_dofs,bc_dofs_x,'rows');
    end
    if (~chk2)
        bc_dofs = union(bc_dofs,bc_dofs_y,'rows');
    end
    if (ndim==3 && ~chk3)
        bc_dofs = union(bc_dofs,bc_dofs_z,'rows');
    end
    
    list1 = find(bc_dofs(:,2)==0);
    bc_dofs(list1,:) = [];
    
    % OUTPUT
    disp_nodes = horzcat(nodes_at_edge1',bc_disp_amount);
    
    % OUTPUT
    disp_dofs = bc_dofs;
    
    %% EDISP
    bc_dofs = union(bc_dofs,fixed_dofs,'rows');
else
    bc_dofs=fixed_dofs;
end

% if (isfield(define_BC,'CDISP'))
%     
%     CDISP=define_BC.CDISP; %=[2.0 -0.5 1 1];'
%     
%     for tt = 1:size(CDISP,1)
%         node_at_specified_coor(tt,1) = find(ismember(data(:,3:3+ndim-1),CDISP(tt,1:ndim),'rows'));
%         
%         if (CDISP(tt,ndim+1));
%             x_dof_cboun(tt,:) = ndim*(node_at_specified_coor(tt,1)-1)+1;
%             x_mag(tt,:)       = CDISP(tt,ndim+1);
%         end;
%         if (CDISP(tt,ndim+2));
%             y_dof_cboun(tt,:) = ndim*(node_at_specified_coor(tt,1)-1)+2;
%             y_mag(tt,:)       = CDISP(tt,ndim+2);
%         end;
%         if (ndim==3 && CDISP(tt,ndim+3));
%             z_dof_cboun(tt,:) = ndim*(node_at_specified_coor(tt,1)-1)+3;
%             z_mag(tt,:)       = CDISP(tt,ndim+3);
%         end;
%     end
%     
%     chk4 = isempty(x_dof_cboun);
%     chk5 = isempty(y_dof_cboun);
%     if (ndim==3)
%         chk6=isempty(z_dof_cboun);
%     end
%     
%     fixed_dofs_cboun = []; % initilize
%     fixed_mag_cboun  = []; % initilize
% 
%     if (~chk4)
%         fixed_dofs_cboun = union(fixed_dofs_cboun,x_dof_cboun');
%         fixed_mag_cboun  = union(fixed_mag_cboun,x_mag');
%     end
%     if (~chk5)
%         fixed_dofs_cboun = union(fixed_dofs_cboun,y_dof_cboun');
%         fixed_mag_cboun  = union(fixed_mag_cboun,y_mag');
%     end
%     if (ndim==3 && ~chk6)
%         fixed_dofs_cboun = union(fixed_dofs_cboun,z_dof_cboun');
%         fixed_mag_cboun  = union(fixed_mag_cboun,z_mag');
%     end
%     
%     fixed_nodes_cboun = node_at_specified_coor;
%     
%     forced_dofs = union(forced_dofs,fixed_dofs_cboun);
%     force_mag   = union(force_mag,fixed_mag_cboun);
%     
% end

% We have bc_dofs and fixed_dofs.
% OUTPUT
% bc_dofs=union(bc_dofs,fixed_dofs,'rows');


end


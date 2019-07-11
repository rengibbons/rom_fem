function [K,R,history_0,history_1] = getKf(sys_vars,ul,old_ul,xl,...
                                           connect_list,element_dof,ndofs,...
                                           sparse_sum,history_0,history_1)

enable_numerical_tangent = sys_vars.enable_numerical_tangent;
physics                  = sys_vars.physics;
n_element_total         = size(connect_list,1);
% n_element_total          = sys_vars.total_element_no; 

K = sparse(ndofs,ndofs);
R = sparse(ndofs,1);

%% LOOP ALL OVER ELEMENTS
% SERIAL
for ie = 1 : n_element_total
    
    % ELEMENT FORMULATION WILL COME HERE.
    % INPUTS: ul: deformation vector
    %       : xl: reference coordinates
    if strcmp(physics,'elastic')
        [K_el,fe] = elmt13(sys_vars,ul(connect_list(ie,4:end)',:),...
                           xl(connect_list(ie,4:end)',3:end));
    elseif  strcmp(physics,'elastic_neoh')
        [K_el,fe] = elmt_elast_neoh(sys_vars,ul(connect_list(ie,4:end)',:),...
                                    xl(connect_list(ie,4:end)',3:end));
    elseif strcmp(physics,'thermo_elastic')
        history_1.x_curr_fixed_el = history_1.x_curr_fixed(connect_list(ie,4:end)',:);
        [K_el,fe,history_1] = elmt_thermo_elast(sys_vars,ul(connect_list(ie,4:end)',:),...
                                                old_ul(connect_list(ie,4:end)',:),...
                                                xl(connect_list(ie,4:end)',3:end),history_0,history_1);
    else
        fprintf('elements physics is unknown. Verify physics in input file.\n')
        stop;
    end
    
    %% ASSEMBLY PART
    [row_map,col_map,index_dofs] = node_mapping_matrix(element_dof(ie,:));
    if enable_numerical_tangent
        K_el_num = numerical_tangent(sys_vars,ul,xl,old_ul,connect_list,ie,history_0,history_1);
        K         = sparse_sum(K,sparse(row_map,col_map,K_el_num,ndofs,ndofs));
        fprintf('\n\nNorm stiff matrix diff = %0.4f\n',norm(K_el-K_el_num)/norm(K_el))
    else
        K = sparse_sum(K,sparse(row_map,col_map,K_el,ndofs,ndofs));
    end
    R = sparse_sum(R,sparse(index_dofs,1,fe,ndofs,1));
end

end

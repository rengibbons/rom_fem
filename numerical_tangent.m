%%
function K_el_num = numerical_tangent(sys_vars,ul,xl,old_ul,connect_list,ie,history_0,history_1)

du = 1e-8;

physics     = sys_vars.physics;
n_dim       = sys_vars.ndim;
n_node_p_el = sys_vars.node_p_el;

if strcmp(physics,'elastic') || strcmp(physics,'elastic_neoh')
    n_dof_p_node = n_dim;
elseif strcmp(physics,'thermo_elastic')
    n_dof_p_node = n_dim+1;
end

n_dof_p_element = n_dof_p_node * n_node_p_el;

K_el_num = zeros(n_dof_p_element);

ul_element = ul(connect_list(ie,4:end)',:);

for ii = 1 : n_dof_p_element

    i_dof = floor((ii-1)/n_dof_p_node) + 1;
    j_dof = mod(ii-1,n_dof_p_node) + 1;

    ul_element(i_dof,j_dof) = ul_element(i_dof,j_dof) + du;
    R1 = get_residual(sys_vars,ul_element,old_ul,xl,history_0,history_1,connect_list,ie);
    
    ul_element(i_dof,j_dof) = ul_element(i_dof,j_dof) - 2*du;
    R2 = get_residual(sys_vars,ul_element,old_ul,xl,history_0,history_1,connect_list,ie);
    
    ul_element(i_dof,j_dof) = ul_element(i_dof,j_dof) + du;
    
    K_el_num(:,ii) = (R2 - R1) / (2 * du);
end

end

%%
function R = get_residual(sys_vars,ul_element,old_ul,xl,history_0,history_1,connect_list,ie)

physics = sys_vars.physics;
if strcmp(physics,'elastic')
    [~,R] = elmt13(sys_vars,ul_element,xl(connect_list(ie,4:end)',3:end));
elseif strcmp(physics,'elastic_neoh')
    [~,R] = elmt_elast_neoh(sys_vars,ul_element,xl(connect_list(ie,4:end)',3:end));
elseif strcmp(physics,'thermo_elastic')
    [~,R,~] = elmt_thermo_elast(sys_vars,ul_element,old_ul(connect_list(ie,4:end)',:),xl(connect_list(ie,4:end)',3:end),history_0(ie,:),history_1(ie,:));
else
    fprintf('elements physics is unknown. Verify physics in input file.\n')
    stop;
end

end
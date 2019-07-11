function[row_map col_map element_connect] = node_mapping_matrix(element_connect)
element_dof_no=length(element_connect);
row_map		= element_connect'*ones(1,element_dof_no);
col_map		= ones(element_dof_no,1)*element_connect;
end
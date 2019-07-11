function [p,ul_gappy,tol_check] = compute_model_3_search_direction(sys_vars,A,B,J,Phi_solution_gappy,connect_list_gappy,element_dof_gappy,ul_gappy0,ul_reduced)

ul_gappy = ul_gappy0 + Phi_solution_gappy * ul_reduced;

[K,R,history_0,history_1] = getKf_model_3(sys_vars,ul,old_ul,xl,...
                                           connect_list,element_dof,ndofs,...
                                           sparse_sum,history_0,history_1,J)

B_bar = B * R_gappy;
A_bar = A * KxPhi_gappy;

[QA,RA] = qr(A_bar);

tol_check = QA' * B_bar;

p = pinv(RA) * QA' * B_bar;
p=0;
ul_gappy=0;
tol_check=0;
end
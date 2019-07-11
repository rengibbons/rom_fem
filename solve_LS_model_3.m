function [ul_reduced,tol_check] = solve_LS_model_3(K_gappy,R_gappy,J,bc_dofs_hdm,...
                                                   ul_reduced,u_increment_bc_dofs_hdm,...
                                                   Phi_solution,A,B,count)

if count==1
    R_gappy = R_gappy - K_gappy(:,bc_dofs_hdm(:,1)) * u_increment_bc_dofs_hdm(:,2);
end

%% SOLVE SYSTEM
K_gappy(:,bc_dofs_hdm(:,1)) = [];
K_gappy(bc_dofs_hdm(:,1),:) = [];
R_gappy(bc_dofs_hdm(:,1))   = [];

R_gappy     = R_gappy(J);
KxPhi_gappy = K_gappy(J,:) * Phi_solution;

B_bar = B * R_gappy;
A_bar = A * KxPhi_gappy;

[QA,RA] = qr(A_bar);

tol_check = QA' * B_bar;

solved_system = pinv(RA) * QA' * B_bar;

ul_reduced = ul_reduced + solved_system;

end
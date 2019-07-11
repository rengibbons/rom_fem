function [A, B] = construct_large_online_matrices(Phi_residual,Phi_stiffness,J)
%% Following Algorithm 2 from Fahrat 2010

[QR,RR] = qr(Phi_residual(J, :),0);
[QK,RK] = qr(Phi_stiffness(J,:),0);

D = pinv(RR) * QR';
E = pinv(RK) * QK';

C = Phi_stiffness * E;
[QC,A] = qr(C,0);

B = QC' * Phi_residual * D;

end
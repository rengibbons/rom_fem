function [cto_mech_voigt,cto_ther_voigt,sigma_voigt,sigma] = mat_thermo_elast_neoh(sys_vars,F,J,T_new,debug,file01)

n_dim = sys_vars.n_dim;
kappa = sys_vars.mat_para.kappa;
G     = sys_vars.mat_para.G;
T0    = sys_vars.mat_para.T0;
alpha = sys_vars.mat_para.alpha;

%% Geometric calculations
C = F'*F;

[N,lambda_sq] = eig(C);
% [N,lambda_sq] = eigs3_fast(C);
lambda_bar    = J^(-1/3)*sqrt(lambda_sq);

first_invariant_bar = lambda_bar(1,1)^2 ...
                    + lambda_bar(2,2)^2 ...
                    + lambda_bar(3,3)^2;

J_th = (1 + alpha*(T_new-T0)) ^ 3;
J_el = J / J_th;

%% Stress (spatial)
A1 = kappa/J_th * (J_el - 1);
A2 = J^(-5/3) * G;

sigma_iso = zeros(3);
n = F*N;
for ii = 1 : 3
    sigma_iso = sigma_iso ...
              + A2*1/lambda_bar(ii,ii)^2 ...
              * (lambda_bar(ii,ii)^2 - 1/3 * first_invariant_bar) ...
              * n(:,ii)*n(:,ii)';
end

sigma = A1*eye(3) + sigma_iso;

%% CTO mechanical push forward
voigt_size               = 3*(n_dim-1);
IxI_voigt                = zeros(voigt_size);
IxI_voigt(1:n_dim,1:n_dim) = ones(n_dim);
II_voigt                 = 1/2*eye(voigt_size);
II_voigt(1:n_dim,1:n_dim)  = eye(n_dim);

% B1 = kappa * (2*J_el - 1) * J_el;
% B2 = 2 * kappa * (J_el - 1) * J_el;
% B3 = 2/3 * J^(-2/3) * G * (C(1,1) + C(2,2) + C(3,3));

% new 
B1 = kappa * (2*J_el - 1) / J_th;
B2 = 2 * kappa * (J_el - 1) / J_th;
B3 = 2/3 * J^(-5/3) * G * (C(1,1) + C(2,2) + C(3,3));

cto_mech_voigt = (B1-1/3*B3)*IxI_voigt + (-B2+B3)*II_voigt;

sigma       = sigma(1:n_dim,1:n_dim);
sigma_voigt = zeros(voigt_size,1);
if n_dim==2
    sigma_voigt(1) = sigma(1,1);
    sigma_voigt(2) = sigma(2,2);
    sigma_voigt(3) = sigma(1,2);
    
%     cto_mech_voigt(1,1) = cto_mech_voigt(1,1) - 4/3*J *  sigma_iso(1,1);
%     cto_mech_voigt(1,2) = cto_mech_voigt(1,2) - 2/3*J * (sigma_iso(1,1)+sigma_iso(2,2));
%     cto_mech_voigt(1,3) = cto_mech_voigt(1,3) - 2/3*J *  sigma_iso(1,2);
%     cto_mech_voigt(2,2) = cto_mech_voigt(2,2) - 4/3*J *  sigma_iso(2,2);
%     cto_mech_voigt(2,3) = cto_mech_voigt(2,3) - 2/3*J *  sigma_iso(1,2);

    % new
    cto_mech_voigt(1,1) = cto_mech_voigt(1,1) - 4/3 *  sigma_iso(1,1);
    cto_mech_voigt(1,2) = cto_mech_voigt(1,2) - 2/3 * (sigma_iso(1,1)+sigma_iso(2,2));
    cto_mech_voigt(1,3) = cto_mech_voigt(1,3) - 2/3 *  sigma_iso(1,2);
    cto_mech_voigt(2,2) = cto_mech_voigt(2,2) - 4/3 *  sigma_iso(2,2);
    cto_mech_voigt(2,3) = cto_mech_voigt(2,3) - 2/3 *  sigma_iso(1,2);
    cto_mech_voigt(2,1) = cto_mech_voigt(1,2);
    cto_mech_voigt(3,1) = cto_mech_voigt(1,3);
    cto_mech_voigt(3,2) = cto_mech_voigt(2,3);
else %n_dim==3
    sigma_voigt(1) = sigma(1,1);
    sigma_voigt(2) = sigma(2,2);
    sigma_voigt(3) = sigma(3,3);
    sigma_voigt(4) = sigma(1,2);
    sigma_voigt(5) = sigma(1,3); 
    sigma_voigt(6) = sigma(2,3);
    
    % TODO: implement cto voigt conversion for 3D
end

%% CTO thermal-mechanical push forward
% cto_ther = 3*alpha*kappa*J_el*J_th^(-1/3)*(1-2*J_el) * eye(voigt_size);

% new
cto_ther = 3*alpha*kappa*J_th^(-4/3)*(1-2*J_el) * eye(voigt_size);

cto_ther_voigt = zeros(voigt_size,1);

if n_dim==2
    cto_ther_voigt(1) = cto_ther(1,1);
    cto_ther_voigt(2) = cto_ther(2,2);
    cto_ther_voigt(3) = cto_ther(1,2);
else %n_dim==3
    cto_ther_voigt(1) = cto_ther(1,1);
    cto_ther_voigt(2) = cto_ther(2,2);
    cto_ther_voigt(3) = cto_ther(3,3);
    cto_ther_voigt(4) = cto_ther(1,2);
    cto_ther_voigt(5) = cto_ther(1,3); 
    cto_ther_voigt(6) = cto_ther(2,3);
end

% change this eventually
% cto_mech_voigt = cto_mech_voigt / J;
% cto_ther_voigt = cto_ther_voigt / J;

end

function [cto_mech_voigt,sigma_voigt,sigma] = mat_elast_neoh(sys_vars,F,J_el,debug,file01)

n_dim = sys_vars.n_dim;
kappa = sys_vars.mat_para.kappa;
G     = sys_vars.mat_para.G;

%% Geometric calculations
if n_dim==2
    F33 = zeros(3);
    F33(1:2,1:2) = F;
    F33(3,3)     = 1;
    F = F33;
end

C = F'*F;

% [N,lambda_sq] = eigs3_fast(C);
[N,lambda_sq] = eig(C);

lambda_bar    = J_el^(-1/3)*sqrt(lambda_sq);

first_invariant_bar = lambda_bar(1,1)^2 ...
                    + lambda_bar(2,2)^2 ...
                    + lambda_bar(3,3)^2;
                
%% Stress (spatial)
A1 = kappa * (J_el - 1);
A2 = J_el^(-5/3) * G;

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

B1 = kappa * (2*J_el - 1);
B2 = 2 * kappa * (J_el - 1);
B3 = 2/3 * J_el^(-5/3) * G * (C(1,1) + C(2,2) + C(3,3));

cto_mech_voigt = (B1-1/3*B3)*IxI_voigt + (-B2+B3)*II_voigt;

sigma       = sigma(1:n_dim,1:n_dim);
sigma_voigt = zeros(1,voigt_size);

if n_dim==2
    sigma_voigt(1) = sigma(1,1);
    sigma_voigt(2) = sigma(2,2);
    sigma_voigt(3) = sigma(1,2);
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

end

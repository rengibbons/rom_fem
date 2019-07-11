function J = determine_greedy_indices(Phi_residual,Phi_stiffness,ni,K_stiffness)
%% Following Algorithm 5 from Fahrat 2010

N     = size(Phi_residual,1);
n_R   = size(Phi_residual,2);
n_K   = size(Phi_stiffness,2);
J     = [];
K_set = [];
n_bar = 0;
m     = 1;

R_vec = Phi_residual(:, 1);
K_vec = Phi_stiffness(:,1);

while n_bar<ni
    k       = setdiff(1:N,J);
    [~,ind] = max(R_vec(k).^2 + K_vec(k).^2);
    J       = union(J,k(ind));

    add_L6 = 0;
    if add_L6
        k = setdiff(1:N,J);
        Gk = g(N,k,K_stiffness);
        GJ = g(N,J,K_stiffness);
        K_set = intersect(Gk,GJ);
        K_set = intersect(K_set,k);
    end
    
    J        = union(J,K_set);
    
    n_bar = n_bar + 1 + size(K_set,1);
    m     = m + 1;
    
    p_R = min([m-1, n_R]);
    p_K = min([m-1, n_K]);
    
    A_R = Phi_residual(J, 1:p_R);
    A_K = Phi_stiffness(J,1:p_K);
    
    phi_Rr = (A_R' * A_R) \ (A_R' * Phi_residual(J, m));
    phi_Kr = (A_K' * A_K) \ (A_K' * Phi_stiffness(J,m));
    
    R_vec = Phi_residual(:, m) - Phi_residual(:, 1:p_R) * phi_Rr;
    K_vec = Phi_stiffness(:,m) - Phi_stiffness(:,1:p_K) * phi_Kr;
    
end

end

%%
function return_set = g(N,J,K_stiffness)
K_red = abs(K_stiffness(J,:));
tol = 1e-12;
[~,return_set] = find(K_red>tol);
return_set = sort(unique(return_set));
end
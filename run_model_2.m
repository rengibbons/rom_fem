function run_model_2(training_folder,ktr,training_params,flags)

fprintf('STARTING MODEL II SIMULATION\n')

% Choose whether to truncate given a specific size
%   0: chooses size of basis on error analysis
%   1: chooses size of basis by specified value below
% choose_trunc_size_soln = 1;

% % Choose an error size
% tol_rom = 1e-12;

load_solution = load(strcat(training_folder,'/model_1_solution_snapshots.mat'));
load_sys_vars = load(strcat(training_folder,'/model_1_sys_vars.mat'));
sys_vars      = load_sys_vars.sys_vars;

sys_vars.model      = 'modelII';
sys_vars.ktr        = ktr;
sys_vars.snap_count = 1;
sys_vars.plot_each_load_step      = flags.plot_each_load_step;
sys_vars.plot_final_configuration = flags.plot_final_configuration;

% Find SVD
[Phi_solution,Sigma_soln,~] = svd(load_solution.solution_snapshots,'econ');

% crude plotting stuff
% format long
% n = size(Sigma_soln,1);
% denom = trace(Sigma_soln.^2);
% err = zeros(n,1);
% sv = zeros(n,1);
% for ii = 1 : n
%     err(ii) = 1 - trace(Sigma_soln(1:ii,1:ii).^2) / denom;
%     sv(ii) = Sigma_soln(1:ii,1:ii);
% end

% scatter(1:n,err)
% [[1:n]' sv err log(err)]
% ktr = 120;
% sys_vars.ktr        = ktr;
% 
% return

% Truncate
sys_vars.Phi_solution = Phi_solution(:,1:ktr);

% podMode  = 1:size(Sigma,1);
% singVals = diag(Sigma);

% nu_training = training_params.nu_training;
% sys_vars.n_training_sims = size(nu_training,1);
% E = sys_vars.mat_para.E;

heat_src_training        = training_params.heat_src_training;
n_training_sims          = size(heat_src_training,1);
sys_vars.n_training_sims = n_training_sims;
n_time_step              = sys_vars.tmax / sys_vars.dt;

for ii = 1 : sys_vars.n_training_sims
    fprintf('\nMODEL II SIMULATION %.0f\n\n',ii)
    sys_vars.i_training_sim = ii;
    
%     % Poission ratio parameter space
%     sys_vars.mat_para.G     = E / (2 * (1 + nu_training(ii)));
%     sys_vars.mat_para.kappa = E / (3 * (1 - 2*nu_training(ii))); 

    % Heat source parameter space
    sys_vars.heat_src.heat_path = zeros(n_time_step,sys_vars.n_dim);
    sys_vars.heat_src.heat_path(:,1) = heat_src_training(ii,1);
    sys_vars.heat_src.heat_path(:,2) = heat_src_training(ii,2);
    
    sys_vars = main_fes_model_2(sys_vars);
end

residual_snapshots  = sys_vars.residual_snapshots(:,1:sys_vars.snap_count-1);
stiffness_snapshots = sys_vars.stiffness_snapshots(:,1:sys_vars.snap_count-1);

% Find SVD
[Phi_residual, Sigma_residual, ~] = svd(residual_snapshots, 'econ');
[Phi_stiffness,Sigma_stiffness,~] = svd(stiffness_snapshots,'econ');

ni = ktr;

J = determine_greedy_indices(Phi_residual,Phi_stiffness,ni,sys_vars.K);

% Truncate
Phi_residual  = Phi_residual(:,1:ktr);
Phi_stiffness = Phi_stiffness(:,1:ktr);

[A, B] = construct_large_online_matrices(Phi_residual,Phi_stiffness,J);

model_2_variables = struct('A',A,...
                           'B',B,...
                           'J',J,...
                           'Phi_solution',sys_vars.Phi_solution,...
                           'ktr',ktr);
                       
save(strcat(training_folder,'/model_2_variables.mat'),'model_2_variables')
% save('gqe_figs/residual_snapshots.mat','residual_snapshots')

end

% % Initialize error array
% error_array = zeros(size(Sigma_soln,1),1);
% 
% % Find where to truncate given an error tolerance, if requested
% denomenator = 0;
% for ii = 1 : size(V_soln,2)
%     denomenator = denomenator + Sigma_soln(ii,ii)^2;
% end
% 
% numerator = 0;
% error_satisfied = false;
% for ii = 1 : size(V_soln,2)
%     numerator = numerator + Sigma_soln(ii,ii)^2;
%     error = 1 - numerator / denomenator;
%     error_array(ii) = error * 100;
%     if error < tol && ~error_satisfied
%         ktr_possible = ii;
%         error_satisfied = true;
%     end
% end
% 
% if ~choose_trunc_size_soln
%     ktr = ktr_possible;
% end

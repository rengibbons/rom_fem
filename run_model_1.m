function run_model_1(input_file,training_folder,training_params,flags)

fprintf('STARTING MODEL I SIMULATION\n')

% Set tolerance
tol = 1.000E-8;

%% OPEN FILES
eval(input_file)
output_folder = 'output';

sys_vars = struct('output_folder',output_folder,...
                  'mesh_input',mesh_input,...
                  'define_BC',define_BC,...
                  'int_rule',int_rule,...
                  'n_dim',n_dim,...
                  'tmax',tmax,...
                  'dt',dt,...
                  'element_type',element_type,...
                  'physics',physics,...
                  'mat_para',mat_para,...
                  'node_p_el',node_p_el,...
                  'total_element_no',0,...
                  'total_node_no',0,...
                  't_para',t_para,...
                  'tol',tol,...
                  'model','modelI',...
                  'debug',flags.debug,...
                  'enable_numerical_tangent',flags.enable_numerical_tangent,...
                  'write_tec',flags.write_tec,...
                  'plot_load',flags.plot_load,...
                  'plot_initial_configuration',flags.plot_initial_configuration,...
                  'plot_each_load_step',flags.plot_each_load_step,...
                  'plot_final_configuration',flags.plot_final_configuration,...
                  'save_solution_vectors',flags.save_solution_vectors);

if strcmp(physics,'thermo_elastic')
    sys_vars.heat_src = heat_src;
end

% nu_training              = training_params.nu_training;
% n_training_sims          = size(nu_training,1);
% sys_vars.n_training_sims = n_training_sims;

heat_src_training        = training_params.heat_src_training;
n_training_sims          = size(heat_src_training,1);
sys_vars.n_training_sims = n_training_sims;
n_time_step              = tmax / dt;

sys_vars.snap_count      = 1;

% Loop over the training parameters
for ii = 1 : n_training_sims
    fprintf('\nMODEL I SIMULATION %.0f\n\n',ii)
    sys_vars.i_training_sim = ii;
    
%     % Poission ratio parameter space
%     sys_vars.mat_para.G     = E / (2 * (1 + nu_training(ii)));
%     sys_vars.mat_para.kappa = E / (3 * (1 - 2*nu_training(ii)));

    % Heat source parameter space
    sys_vars.heat_src.heat_path = zeros(n_time_step,n_dim);
    sys_vars.heat_src.heat_path(:,1) = heat_src_training(ii,1);
    sys_vars.heat_src.heat_path(:,2) = heat_src_training(ii,2);
    
    sys_vars = main_fes_model_1(sys_vars);
end

% Save solution vectors
solution_snapshots = sys_vars.solution_snapshots;

warning('off','all')
mkdir(training_folder)
warning('on','all')

save(strcat(training_folder,'/model_1_solution_snapshots.mat'),'solution_snapshots')

% Clear solution vector from sys_vars to make structure smaller.
sys_vars.solution_snapshots = 0;

save(strcat(training_folder,'/model_1_sys_vars.mat'),'sys_vars')
end

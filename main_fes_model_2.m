function sys_vars = main_fes_model_2(sys_vars)

mesh_input               = sys_vars.mesh_input;
define_BC                = sys_vars.define_BC;
dt                       = sys_vars.t_para.dt;
tmax                     = sys_vars.t_para.tmax;
physics                  = sys_vars.physics;
n_dim                    = sys_vars.n_dim;
tol_fem                  = sys_vars.tol;
plot_each_load_step      = sys_vars.plot_each_load_step;
plot_final_configuration = sys_vars.plot_final_configuration;
Phi_solution             = sys_vars.Phi_solution;
ktr                      = sys_vars.ktr;

%% MESH GEOMETRY
[sys_vars,xl,connect_list,element_dof,dof_list] = mesh(sys_vars,mesh_input);

%% SET BOUNDARIES - GET FIXED DEGREE OF FREEDOMS
[fixed_dofs,fixed_nodes] = set_dirichlet_boundary(sys_vars,mesh_input,define_BC,dof_list);

%% SET DISPLACEMENT LOADING (force loading not implememnted currently)
[bc_dofs_hdm,~,disp_nodes] = set_disp_boundary(sys_vars,mesh_input,define_BC,fixed_dofs);
[force_dofs,force_mag] = set_force_boundary(sys_vars,mesh_input,define_BC);

%% INITIALIZE EVERYTHING
history_0 = struct();
history_1 = struct();
if strcmp(physics,'elastic') || strcmp(physics,'elastic_neoh')
    ndofs         = sys_vars.total_node_no * sys_vars.n_dim;
    ndof_per_node = sys_vars.n_dim;
elseif strcmp(physics,'thermo_elastic')
    ndofs         = sys_vars.total_node_no * (sys_vars.n_dim + 1);
    ndof_per_node = sys_vars.n_dim + 1;
end
sparse_sum = @plus;

x_curr_hdm = xl(:,3:3+ndof_per_node-1);
x_reff_hdm = xl(:,3:3+ndof_per_node-1);
if strcmp(physics,'thermo_elastic')
    x_reff_hdm(:,end) = sys_vars.mat_para.T0;
    x_curr_hdm(:,end) = sys_vars.mat_para.T0;
end
x_tec             = horzcat(x_reff_hdm,x_curr_hdm);
u_inc_bc_dofs_hdm = bc_dofs_hdm;
ul_hdm            = zeros(sys_vars.total_node_no,ndof_per_node);
ul_list_hdm       = zeros(ndofs,1);
ul_list_red   = zeros(ktr,1);
old_ul_hdm        = zeros(sys_vars.total_node_no,ndof_per_node);
old_ul_list_hdm   = zeros(ndofs,1);
% old_ul_list_red   = zeros(ktr,1);
free_dofs_hdm     = setdiff(1:ndofs,bc_dofs_hdm(:,1))';
n_free_dofs       = size(free_dofs_hdm,1);

%% LOADING FUNCTION DEFINITION
% MONOTONIC LOADING t_0=[0,0] t_1=[1,1] t_max=[max,1]
n_time_step     = tmax / dt;
sys_vars.n_time_step = n_time_step;
simulation_time = 0 : dt : tmax;

if sys_vars.i_training_sim==1
    n_iter_guess = 10;
    sys_vars.residual_snapshots  = zeros(n_free_dofs,n_iter_guess*n_time_step*sys_vars.n_training_sims);
    sys_vars.stiffness_snapshots = zeros(n_free_dofs,n_iter_guess*n_time_step*sys_vars.n_training_sims);
end

%% FIRST LOOP: START -> LOAD/TIME STEP
global_timer = tic;
for it = 1 : n_time_step
    fprintf('STEP: %5d    sim time: %6.3f\n',it,simulation_time(1+it));
    
%     sys_vars.time_step = it;
    
    load_factor = simulation_time(1+it)/tmax;
    if load_factor>1.-eps()
        load_factor = 1;
    end
    
    %% COMPUTE INCREMENTAL LOADING
    u_inc_bc_dofs_hdm(:,2) = bc_dofs_hdm(:,2) / n_time_step; 
    
    % SET INITIAL RESIDUAL
    count = 0;
    
    if strcmp(physics,'thermo_elastic')
        history_1.x_curr_fixed = x_curr_hdm;
        history_1.it           = it;
    end
    
    %% GLOBAL NEWTON ITERATION
    while (1)
        count = count+1;
        
        [K,R,history_0,history_1] = getKf(sys_vars,ul_hdm,old_ul_hdm,xl,...
                                          connect_list,element_dof,ndofs,...
                                          sparse_sum,history_0,history_1);

        %% SOLVE THE SYSTEM
        R(force_dofs) = R(force_dofs) + load_factor.*force_mag;
        
        if count==1
            R = R - K(:,bc_dofs_hdm(:,1)) * u_inc_bc_dofs_hdm(:,2);
        end

        K(:,bc_dofs_hdm(:,1)) = [];
        K(bc_dofs_hdm(:,1),:) = [];
        R(bc_dofs_hdm(:,1))   = [];
        
        KxPhi         = K * Phi_solution;
        solved_system = (KxPhi' * KxPhi) \ (KxPhi' * R);
        
        
%         ul_list_red                     = old_ul_list_red + solved_system;
%         ul_list_hdm(free_dofs_hdm)      = Phi_solution * ul_list_red + ul_list_hdm(free_dofs_hdm); % trying this one
        ul_list_red                = ul_list_red + solved_system; % trying this one 
        ul_list_hdm(free_dofs_hdm) = Phi_solution * solved_system + ul_list_hdm(free_dofs_hdm);

        
        ul_list_hdm(bc_dofs_hdm(:,1)) = u_inc_bc_dofs_hdm(:,2) + old_ul_list_hdm(bc_dofs_hdm(:,1));
        residual                        = norm(R);
        
        if count==1
            residual_initial = residual;
        end

        %% STORE SOLUTION
        ul_hdm(:,1) = ul_list_hdm(1:ndof_per_node:end,1);
        ul_hdm(:,2) = ul_list_hdm(2:ndof_per_node:end,1);
        if (strcmp(physics,'elastic')  || strcmp(physics,'elastic_neoh')) && n_dim == 3 
            ul_hdm(:,3) = ul_list_hdm(3:ndof_per_node:end,1);
        elseif strcmp(physics,'thermo_elastic')
            ul_hdm(:,3) = ul_list_hdm(3:ndof_per_node:end,1);
            if n_dim==3
                ul_hdm(:,4) = ul_list_hdm(4:ndof_per_node:end,1);
            end
        end
        
        if sys_vars.snap_count==1
            sys_vars.K = K;
        end

        % Save snapshots
        model = 1;
        sys_vars.residual_snapshots(:,sys_vars.snap_count)  = R;
        if model==1
            sys_vars.stiffness_snapshots(:,sys_vars.snap_count) = R;
        elseif model==2
            sys_vars.stiffness_snapshots(:,sys_vars.snap_count) = KxPhi * solved_system;
        end
        sys_vars.snap_count = sys_vars.snap_count + 1;

        %% OUTPUT TO FILE RESIDUAL
        fprintf('Iter: %3d \t r: %5.2e \t N(r): %5.2e\n',count,residual,residual/residual_initial);

        %% EXIT LOOP IF SATISFIED (NUMBER OF STEPS CAN BE ADDED TOO.)
        if abs(residual)<tol_fem
            old_ul_hdm      = ul_hdm;
            old_ul_list_hdm = ul_list_hdm;
%             old_ul_list_red = ul_list_red;
            x_curr_hdm      = x_reff_hdm + ul_hdm;
            x_tec(:,ndof_per_node+1:2*ndof_per_node) = x_curr_hdm;
            history_0       = history_1;
            fprintf('\n');
            break;
        elseif count==50
            fprintf('\n NO CONVERGENCE, residual= %6.4f \n \n',residual);
            stop;
        end
    end

    %% PLOT AT EACH CONVERGED SOLUTION FOR EACH LOAD STEP.
    if plot_each_load_step
        plot_mesh(sys_vars,x_curr_hdm,connect_list,fixed_nodes(:,1),disp_nodes)
    end
end
toc(global_timer)

%% PLOT FIGURE;
if plot_final_configuration
    plot_mesh(sys_vars,x_curr_hdm,connect_list,fixed_nodes(:,1),disp_nodes)
end

end

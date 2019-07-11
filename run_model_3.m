function run_model_3(training_folder,flags)

fprintf('STARTING MODEL III SIMULATION\n')

%% LOAD FILES AND VARIABLES
load_sys_vars     = load(strcat(training_folder,'/model_1_sys_vars.mat'));
load_model_2      = load(strcat(training_folder,'/model_2_variables.mat'));
sys_vars          = load_sys_vars.sys_vars;
model_2_variables = load_model_2.model_2_variables;

%% System parameters
sys_vars.model           = 'modelIII';
n_dim                    = sys_vars.n_dim;
mesh_input               = sys_vars.mesh_input;
define_BC                = sys_vars.define_BC;
dt                       = sys_vars.t_para.dt;
tmax                     = sys_vars.t_para.tmax;
physics                  = sys_vars.physics;
n_time_step              = sys_vars.n_time_step;
plot_each_load_step      = flags.plot_each_load_step;
plot_final_configuration = flags.plot_final_configuration;

%% Heat path parameter choice
x0 = 0.1*[0.1,0.1];
xf = 0.9*[0.1,0.1];

x_path    = linspace(x0(1),xf(1),n_time_step);
y_path    = linspace(x0(2),xf(2),n_time_step);
sys_vars.heat_src.heat_path = [x_path' y_path'];


% Model II variables
A            = model_2_variables.A;
B            = model_2_variables.B;
J            = model_2_variables.J;
Phi_solution = model_2_variables.Phi_solution;
ktr          = model_2_variables.ktr;

% % Convert material parameters for nu variation
% nu = 0.44;
% E = sys_vars.mat_para.E;
% sys_vars.mat_para.G     = E / (2 * (1 + nu));
% sys_vars.mat_para.kappa = E / (3 * (1 - 2*nu));

%% MESH GEOMETRY
[sys_vars,xl_hdm,connect_list,element_dof,dof_list] = mesh(sys_vars,mesh_input);

%% SET BOUNDARIES - GET FIXED DEGREE OF FREEDOMS
[fixed_dofs_hdm,fixed_nodes_hdm] = set_dirichlet_boundary(sys_vars,mesh_input,define_BC,dof_list);

%% SET DISPLACEMENT LOADING (force loading not implememnted currently)
[bc_dofs_hdm,~,disp_nodes_hdm] = set_disp_boundary(sys_vars,mesh_input,define_BC,fixed_dofs_hdm);
[force_dofs_hdm,force_mag_hdm] = set_force_boundary(sys_vars,mesh_input,define_BC);

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

x_curr_hdm = xl_hdm(:,3:3+ndof_per_node-1);
x_reff_hdm = xl_hdm(:,3:3+ndof_per_node-1);
if strcmp(physics,'thermo_elastic')
    x_reff_hdm(:,end)   = sys_vars.mat_para.T0;
end
x_tec             = horzcat(x_reff_hdm,x_curr_hdm);
u_inc_bc_dofs_hdm = bc_dofs_hdm;
ul_red            = zeros(ktr,1);
ul_hdm            = zeros(sys_vars.total_node_no,ndof_per_node);
ul_list_hdm       = zeros(ndofs,1);
old_ul_hdm        = zeros(sys_vars.total_node_no,ndof_per_node);
old_ul_list_hdm   = zeros(ndofs,1);
free_dofs_hdm     = setdiff(1:ndofs,bc_dofs_hdm(:,1))';
sys_vars.total_element_no = size(element_dof,1);

% Find the gappy elements needed for assembly.
J_dof_hdm      = free_dofs_hdm(J);
elements_gappy = [];
for ii = 1 : size(J_dof_hdm,1)
    [elements_gappy_new,~] = find(element_dof==J_dof_hdm(ii));
    elements_gappy = union(elements_gappy,elements_gappy_new);
end
connect_list_gappy = connect_list(elements_gappy,:);
element_dof_gappy  = element_dof(elements_gappy,:);

%% LOADING FUNCTION DEFINITION
simulation_time = 0 : dt : tmax;

tol = 1e-8;

%% LOOP OVER TIME STEPS
tic
for it = 1 : n_time_step
    fprintf('\nSTEP: %5d    sim time: %6.3f\n',it,simulation_time(1+it));
    
    count = 0;
    
    if strcmp(physics,'thermo_elastic')
        history_1.x_curr_fixed = x_curr_hdm;
        history_1.it           = it;
    end

    load_factor = simulation_time(1+it)/tmax;
    if load_factor>1.-eps()
        load_factor = 1;
    end
    
    %% COMPUTE INCREMENTAL LOADING
    u_inc_bc_dofs_hdm(:,2) = bc_dofs_hdm(:,2) / n_time_step;
    
    %% GLOBAL NEWTON ITERATION
    while (1)
        count = count + 1;

        [K_gappy,R_gappy,history_0,history_1] = getKf(sys_vars,ul_hdm,old_ul_hdm,xl_hdm,...
                                                      connect_list_gappy,element_dof_gappy,ndofs,...
                                                      sparse_sum,history_0,history_1);

        R_gappy(force_dofs_hdm) = R_gappy(force_dofs_hdm) + load_factor * force_mag_hdm;
        
        %% SOLVE SYSTEM
        [ul_red,tol_check] = solve_LS_model_3(K_gappy,R_gappy,J,bc_dofs_hdm,...
                                                  ul_red,u_inc_bc_dofs_hdm,...
                                                  Phi_solution,A,B,count);
      
        ul_list_hdm(free_dofs_hdm)      = Phi_solution * ul_red;
        ul_list_hdm(bc_dofs_hdm(:,1),:) = u_inc_bc_dofs_hdm(:,2) + old_ul_list_hdm(bc_dofs_hdm(:,1),1);
      
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
        
        fprintf('Iter: %3d \t QTB: %5.2e \t r %5.2e\n',count,abs(norm(tol_check)),norm(R_gappy));
        
        %% EXIT LOOP IF SATISFIED
        if abs(norm(tol_check))<tol
            old_ul_red = ul_red;
            old_ul_hdm       = ul_hdm;
            old_ul_list_hdm  = ul_list_hdm;
            x_curr_hdm       = x_reff_hdm + ul_hdm;
            history_0    = history_1;
            break;
        elseif count==50
            fprintf('\n NO CONVERGENCE, residual= %6.4f \n \n',residual);
            error('NO CONVERGENCE')
        end      
    end 
    
    %% PLOT AT EACH CONVERGED SOLUTION FOR EACH LOAD STEP.
    if plot_each_load_step
        plot_mesh(sys_vars,x_curr_hdm,connect_list,fixed_nodes_hdm(:,1),disp_nodes_hdm)
    end
end
toc

%% PLOT FIGURE;
if plot_final_configuration
    plot_mesh(sys_vars,x_curr_hdm,connect_list,fixed_nodes_hdm(:,1),disp_nodes_hdm)
    print('gqe_figs/therel_ht_path_ktr150','-depsc','-tiff')
    save(strcat('gqe_figs/ul_ktr150.mat'),'ul_hdm')
end

% dlmwrite('elast_comp/m3_044',ul_list_hdm,'delimiter','\t','precision',10)
    
end



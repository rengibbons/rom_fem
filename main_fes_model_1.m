function sys_vars = main_fes_model_1(sys_vars)

%% GLOBAL VARIABLES
output_folder              = sys_vars.output_folder;
mesh_input                 = sys_vars.mesh_input;
define_BC                  = sys_vars.define_BC;
dt                         = sys_vars.t_para.dt;
tmax                       = sys_vars.t_para.tmax;
physics                    = sys_vars.physics;
n_dim                      = sys_vars.n_dim;
tol                        = sys_vars.tol;
write_tec                  = sys_vars.write_tec;
plot_load                  = sys_vars.plot_load;
plot_initial_configuration = sys_vars.plot_initial_configuration;
plot_each_load_step        = sys_vars.plot_each_load_step;
plot_final_configuration   = sys_vars.plot_final_configuration;
save_solution_vectors      = sys_vars.save_solution_vectors;

output_file         = strcat(output_folder,'/output.dat');
output_tecplot_file = strcat(output_folder,'/output_tecplot.dat');
file01              = fopen(output_file,'w');
file02              = fopen(output_tecplot_file,'w');

sys_vars.file01 = file01;
sys_vars.file02 = file02;

%% WRITE TITLE TO FILES
fprintf(file01,'\n FES - Finite Element Solver \n ');
fprintf(file01,'MATLAB version 3/23/2015 \n ');
fprintf(file01,'v2 mixed brick 8 element is added: 04/25/2015 \n ');
fprintf(file01,'v2 bug how to calculate deformation gradient: 04/21/2015 \n ');
fprintf(file01,'v2 1 point integration is added: 04/20/2015 \n ');

%% MESH GEOMETRY
[sys_vars,xl,connect_list,element_dof,dof_list] = mesh(sys_vars,mesh_input);

%% SET BOUNDARIES - GET FIXED DEGREE OF FREEDOMS
[fixed_dofs,fixed_nodes] = set_dirichlet_boundary(sys_vars,mesh_input,define_BC,dof_list);

%% SET DISPLACEMENT LOADING (force loading not implememnted currently)
[bc_dofs,~,disp_nodes] = set_disp_boundary(sys_vars,mesh_input,define_BC,fixed_dofs);
[force_dofs,force_mag] = set_force_boundary(sys_vars,mesh_input,define_BC);

%% PLOT FIGURE;
if (plot_initial_configuration)
    plot_mesh(sys_vars,xl,connect_list,fixed_nodes(:,1),disp_nodes)
end

%% FOR TECPLOT
if (write_tec) 
    fprintf(file02,'TITLE = "  "\n');
    if (n_dim==3)
        fprintf(file02,'VARIABLES = "x_mat","y_mat","z_mat","x_cur","y_cur","z_cur"\n');
    elseif (n_dim==2)
        fprintf(file02,'VARIABLES = "x_mat","y_mat","x_cur","y_cur"\n');
    end
end

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

x_curr = xl(:,3:3+ndof_per_node-1);
x_reff = xl(:,3:3+ndof_per_node-1);
if strcmp(physics,'thermo_elastic')
    x_reff(:,end)   = sys_vars.mat_para.T0;
    x_curr(:,end)   = sys_vars.mat_para.T0;
end
x_tec         = horzcat(x_reff,x_curr);
u_inc_bc_dofs = bc_dofs;
ul            = zeros(sys_vars.total_node_no,ndof_per_node);
ul_list       = zeros(ndofs,1);
old_ul        = zeros(sys_vars.total_node_no,ndof_per_node);
old_ul_list   = zeros(ndofs,1);
free_dofs     = setdiff(1:ndofs,bc_dofs(:,1))';
n_free_dofs   = size(free_dofs,1);

%% LOADING FUNCTION DEFINITION
% MONOTONIC LOADING t_0=[0,0] t_1=[1,1] t_max=[max,1]
n_time_step          = tmax / dt;
sys_vars.n_time_step = n_time_step;
simulation_time      = 0 : dt : tmax;

if plot_load
    plot_loading_trend([0 0],tmax,n_time_step,n_time_step)
end

if save_solution_vectors && sys_vars.i_training_sim==1
    sys_vars.solution_snapshots = zeros(n_free_dofs,n_time_step*sys_vars.n_training_sims);
end

%% FIRST LOOP: START -> LOAD/TIME STEP
global_timer = tic;
for it = 1 : n_time_step
    fprintf(file01,'\n STEP: %5d',it);
    fprintf('STEP: %5d    sim time: %6.3f\n',it,simulation_time(1+it));

    load_factor = simulation_time(1+it)/tmax;
    if load_factor>1.-eps()
        load_factor = 1;
    end

    % TECPLOT
    if (write_tec)
        if (n_dim==3)
            fprintf(file02,'ZONE N=%d, E= %d, F=FEPOINT, ET=BRICK, SOLUTIONTIME= %8.4e\n',sys_vars.total_node_no,sys_vars.total_element_no,simulation_time(1+it));
            fprintf(file02,'%10.4e %10.4e %10.4e %10.4e %10.4e %10.4e\n',x_tec');
            fprintf(file02,'%g %g %g %g %g %g %g %g\n',connect_list(:,4:4+node_p_el-1)');
        elseif (n_dim==2)
            fprintf(file02,'ZONE N=%d, E= %d, F=FEPOINT, ET=QUADRILATERAL, SOLUTIONTIME= %8.4e\n',sys_vars.total_node_no,sys_vars.total_element_no,simulation_time(1+it));
            fprintf(file02,'%10.4e %10.4e %10.4e %10.4e\n',x_tec');
            fprintf(file02,'%g %g %g %g\n',connect_list(:,4:4+sys_vars.node_p_el-1)');
        end
    end
    
    %% COMPUTE INCREMENTAL LOADING
    u_inc_bc_dofs(:,2) = bc_dofs(:,2) / n_time_step;
    
    % SET INITIAL RESIDUAL
    count = 0;
    
    if strcmp(physics,'thermo_elastic')
        history_1.x_curr_fixed = x_curr;
        history_1.it           = it;
    end
    
    %% GLOBAL NEWTON ITERATION
    while (1)
        count = count+1;
        
        [K,R,history_0,history_1] = getKf(sys_vars,ul,old_ul,xl,...
                                          connect_list,element_dof,ndofs,...
                                          sparse_sum,history_0,history_1);

        %% SOLVE THE SYSTEM
        R(force_dofs) = R(force_dofs) + load_factor.*force_mag;
        [ul_list,R] = solve_LS(ul_list,old_ul_list,K,R,u_inc_bc_dofs,ndofs,count);
        residual = compute_residual(R);
        if count==1
            residual_initial = residual;
        end

        ul(:,1) = ul_list(1:ndof_per_node:end,1);
        ul(:,2) = ul_list(2:ndof_per_node:end,1);
        if (strcmp(physics,'elastic')  || strcmp(physics,'elastic_neoh')) && n_dim == 3 
            ul(:,3) = ul_list(3:ndof_per_node:end,1);
        elseif strcmp(physics,'thermo_elastic')
            ul(:,3) = ul_list(3:ndof_per_node:end,1);
            if n_dim==3
                ul(:,4) = ul_list(4:ndof_per_node:end,1);
            end
        end
        
        %% OUTPUT TO FILE RESIDUAL
        fprintf(file01,'\nIter: %3d \t r: %5.2e \t N(r): %5.2e \t time: %5.2e ',count,residual,residual/residual_initial,toc(global_timer));
        fprintf('Iter: %3d \t r: %5.2e \t N(r): %5.2e \t time: %5.2e\n',count,residual,residual/residual_initial,toc(global_timer));

        %% EXIT LOOP IF SATISFIED
        if abs(residual)<tol
            old_ul       = ul;
            old_ul_list  = ul_list;
            x_curr       = x_reff + ul;
            x_tec(:,ndof_per_node+1:2*ndof_per_node) = x_curr;
            history_0    = history_1;
            fprintf('\n');
            if save_solution_vectors
                % Save solution for training an HDM.
                sys_vars.solution_snapshots(:,sys_vars.snap_count) = ul_list(free_dofs);
                sys_vars.snap_count = sys_vars.snap_count + 1;
            end
            break;
        elseif count==50
            fprintf('\n NO CONVERGENCE, residual= %6.4f \n \n',residual);
            fprintf(file01,'NO CONVERGENCE, residual= %6.4f \n \n',residual);
            error('NO CONVERGENCE')
        end
    end
    
    %% PLOT AT EACH CONVERGED SOLUTION FOR EACH LOAD STEP.
    if plot_each_load_step
        plot_mesh(sys_vars,x_curr,connect_list,fixed_nodes(:,1),disp_nodes)
    end
end
toc(global_timer)

%% PLOT FIGURE;
if plot_final_configuration
    plot_mesh(sys_vars,x_curr,connect_list,fixed_nodes(:,1),disp_nodes)
%     print('gqe_figs/therel_ht_point','-depsc','-tiff')
%     save('gqe_figs/ul_hdm.mat','ul')
end

% dlmwrite('elast_comp/m1_044',full(ul_list),'delimiter','\t','precision',10)

%% CLOSE THE FILE
fclose(file01);
fclose(file02);

end

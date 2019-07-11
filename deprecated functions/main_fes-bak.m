clear all; close all; clc;

%% GLOBAL VARIABLES
rand(5,5)*inv(rand(5,5));

tic
global debug file01 node_p_el int_rule total_element_no ndim total_node_no element_type node_p_elm
%% MOST GLOBAL DESICIONS
mixed =0;
serial =1; % or 0 -> parallel
debug=0;
%% OPEN FILES FOR POST PROSESSING
global_timer=tic;
Output_file='output/output.dat';
Otecplo_file='output/outp_tecplo.dat';
file01=fopen(Output_file,'w');
file02=fopen(Otecplo_file,'w');
%% WRITE TITLE TO FILES

fprintf(file01,'\n FES - Finite Element Solver \n ');
fprintf(file01,'MATLAB version 3/23/2015 \n ');
fprintf(file01,'v2 mixed brick 8 element is added: 04/25/2015 \n ');
fprintf(file01,'v2 bug how to calculate deformation gradient: 04/21/2015 \n ');
fprintf(file01,'v2 1 point integration is added: 04/20/2015 \n ');

%% FORMAT DISPLAY

format shorte
%% FLAGS

write_tec=1; plt_load=0; plt_each_load_step=0; plt_final_configuration=0; plot_initial_configuration=0;
%% INPUT FILE
% INPUT THE COMMON NAME FOR MESH FILES
% SET MATERIAL PARAMETERS

% input_file_brick8_single;
% input_file_brick8;
input_file_quad4;
%  input_file_quad4_single;
%% MESH GEOMETRY

[xl, connect_list,element_dof,dof_list]= Mesh_main(mesh_input); % xl: reference coordinates
%% SET BOUNDARIES - GET FIXDE DEGREE OF FREEDOMS

[fixed_dofs,fixed_nodes]=set_dirichlet_boundary(mesh_input,define_BC,dof_list);
% [fixed_dofs,fixed_nodes]=dirichlet_boundary_set(mesh_input,reference_geometry,dof_list);
%% SET DISPLACEMENT LOADING (force loading not implememnted currently)

[bc_dofs,disp_dofs,disp_nodes]=set_disp_boundary(mesh_input,define_BC,dof_list,fixed_dofs);
[forced_dofs,force_mag]=set_forc_boundary(mesh_input,define_BC,dof_list,fixed_dofs);

% OLD -> [bc_dofs,disp_dofs,disp_nodes]=disp_boundary_set(mesh_input,reference_geometry,dof_list,fixed_dofs);
%% PLOT FIGURE;

if (plot_initial_configuration); PlotMesh(xl,connect_list,fixed_nodes,disp_nodes); end;
%% FOR TECPLOT

if (write_tec) fprintf(file02,'TITLE = "  "\n'); end;
if (write_tec)
    if (ndim==3)
        fprintf(file02,'VARIABLES = "x_mat","y_mat","z_mat","x_cur","y_cur","z_cur"\n');
    elseif (ndim==2)
        fprintf(file02,'VARIABLES = "x_mat","y_mat","x_cur","y_cur"\n');
    end
end

%% INITILIZE EVERYTHING, YO!
n_element=size(connect_list,1);
x_curr=xl(:,3:3+ndim-1);
x_reff=xl(:,3:3+ndim-1); % set material configuration
x_tec=horzcat(x_reff,x_curr);
u_increment_bc_dofs=bc_dofs; % WE NEED THIS SO WE CAN UPDATE u_increment_bc_dofs WITH FACTOR*BC_DOFS
ul=zeros(total_node_no,ndim); % Displacement vector
ul_list=zeros(total_node_no*ndim,1);
old_ul_list=ul_list;
old_dof_list=dof_list;

%% PREALLOCATE/INITLIZE FOR MIXED ELEMENT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mixed=0;
if (mixed)
    i_J = zeros(1,1)+1;
    i_p = zeros(1,1);
    i_inv_k_pj = zeros(1,1);
    i_k_pu = zeros(node_p_elm*ndim,1);
    i_k_JJ = zeros(1,1);
    i_F_p = zeros(1,1);
    i_F_j = zeros(1,1);
    
    i_locate=vertcat(i_J,i_p,i_inv_k_pj,i_k_pu,i_k_JJ,i_F_p,i_F_j);
    add_field_e=zeros(length(i_locate),length(i_locate));
    add_field=repmat(i_locate,1,n_element);
    du=zeros(total_node_no,ndim);
else
    add_field_e=0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

ndofs = total_node_no*ndim;
sparse_sum = @plus; % sum sparse matrices.
K = sparse(total_node_no*ndim,total_node_no*ndim);
R = sparse(total_node_no*ndim,1);
% dof_list;
dof_listn=dof_list;
%% LOADING FUNCTION DEFINITION
% MONOTONIC LOADING t_0=[0,0] t_1=[1,1] t_max=[max,1]
tmax=5;
dt=1;
total_load_step=tmax/dt;
num_load_step=5; %steps
simulation_time=0:dt:tmax;

if (plt_load) plt_loading_trend([0 0],tmax,total_load_step,num_load_step); end;
%% FIRST LOOP: START -> LOAD/TIME STEP

res = zeros(ndofs-size(bc_dofs,1),1);
Fef = zeros(ndofs,1);
Fef(forced_dofs) = force_mag;
dlam = 0.;
sMu=1;
lambda=0;

for it=1:total_load_step
    fprintf(file01,'\n\n STEP: %5d',it);
    d(4,1)=simulation_time(1+it);
    
    load_factor = simulation_time(1+it)/tmax;
    if load_factor>1.-eps()
        load_factor=1;
    end
    
    delul_list  =zeros(ndofs-size(bc_dofs,1),1);
    
    % TECPLOT
    if (write_tec)
        if (ndim==3)
            fprintf(file02,'ZONE N=%d, E= %d, F=FEPOINT, ET=BRICK, SOLUTIONTIME= %8.4e\n',total_node_no,total_element_no,d(4,1));
            fprintf(file02,'%10.4e %10.4e %10.4e %10.4e %10.4e %10.4e\n',x_tec');
            fprintf(file02,'%g %g %g %g %g %g %g %g\n',connect_list(:,4:4+node_p_elm-1)');
        elseif (ndim==2)
            fprintf(file02,'ZONE N=%d, E= %d, F=FEPOINT, ET=QUADRILATERAL, SOLUTIONTIME= %8.4e\n',total_node_no,total_element_no,d(4,1));
            fprintf(file02,'%10.4e %10.4e %10.4e %10.4e\n',x_tec');
            fprintf(file02,'%g %g %g %g\n',connect_list(:,4:4+node_p_elm-1)');
        end
    end
    
    %% COMPUTE INCREMENTAL LOADING
    if ( it<=num_load_step )
        u_increment_bc_dofs(:,2) = dt.*bc_dofs(:,2)/num_load_step;   % scale dirichlet bc's with dt
    else
        % Reaching the maximum iteration step, no increments of loading.
        u_increment_bc_dofs(:,2) = 0;
    end
    % SET TOLERANCE AND INITIAL RESIDUAL
    tol=1.000E-8;
    residual=1;
    residual_0=1;
    count=0;
    %% GLOBAL NEWTON ITERATION
    while (1)
        count=count+1;
        
        K = sparse(total_node_no*ndim,total_node_no*ndim);
        R = sparse(total_node_no*ndim,1);
        %% LOOP ALL OVER ELEMENTS
        % SERIAL
        if (serial)
            for ie=1:total_element_no
                
                % ELEMENT FORMULATION WILL COME HERE.
                % INPUTS: ul: deformation vector
                %       : xl: reference coordinates
                
                [Ke,fe]=elmt13(d,ul(connect_list(ie,4:end)',1:ndim),xl(connect_list(ie,4:end)',3:3+ndim-1));
                %               [Ke,fe,add_field_e]=mixed_elmt13(d,ul(connect_list(ie,4:end)',1:ndim),du(connect_list(ie,4:end)',1:ndim),xl(connect_list(ie,4:end)',3:3+ndim-1),add_field(:,ie));
                %                 [Ke,fe,add_field_e]=mixed_elmt13_no_growth(d,ul(connect_list(ie,4:end)',1:ndim),du(connect_list(ie,4:end)',1:ndim),xl(connect_list(ie,4:end)',3:3+ndim-1),add_field(:,ie));
                
                add_field(:,ie)=add_field_e;
                %% ASSEMBLY PART
                [row_map col_map index_dofs] = node_mapping_matrix(element_dof(ie,:));
                K = sparse_sum(K,sparse(row_map,col_map,Ke,ndofs,ndofs));
                R = sparse_sum(R,sparse(index_dofs,1,fe,ndofs,1));
            end
            
        else
            %% LOOP ALL OVER ELEMENTS (NEEDS IMPROVEMENT)
            % PARALLEL
            parfor ie=1:total_element_no
                %%% ELEMENT FORMULATION WILL COME HERE.
                [Ke,fe]=elmt13_par(d,ul(connect_list(ie,4:end)',1:ndim),xl(connect_list(ie,4:end)',3:3+ndim-1),node_p_el,int_rule,total_node_no,ndim,element_type,node_p_elm,debug,file01);
                [row_map col_map index_dofs] = node_mapping_matrix(element_dof(ie,:));
                
                K = sparse_sum(K,sparse(row_map,col_map,Ke,ndofs,ndofs));
                R = sparse_sum(R,sparse(index_dofs,1,fe,ndofs,1));
                
            end
            %             delete(gcp)
        end
        %% external force
        
        
        %% SOLVE THE SYSTEM
        % Global stiffness K and force F vector are assembled.
        % Now, we need to solve the system -> K d = F for d: incremental deformation
        %
        if (0)
            R(forced_dofs) = R(forced_dofs) + load_factor.*force_mag;
            [ul_list,old_ul_list,R] = solve_LS(ul_list,old_ul_list,K,R,u_increment_bc_dofs,ndofs,count);
            [residual]=compute_residual2(R); if (count==1); residual_initial=residual; end;
        else
            [ul_list,old_ul_list,res,sMu,dlam,delul_list] = solve_arclength(delul_list,ul_list,old_ul_list,K,R,u_increment_bc_dofs,ndofs,count,res,Fef,lambda,dlam,sMu);
            
            
            
            [residual]=compute_residual2(res); if (count==1) residual_initial=residual; end;
        end
        %         [residual]=compute_residual(R,u_increment_bc_dofs(:,1)); if (count==1) residual_initial=residual; end;
        res_list(count,it+1)=residual;
        %% OUTPUT TO FILE RESIDUAL
        
        fprintf(file01,'\nIter: %3d \t r: %5.2e \t N(r): %5.2e \t time: %5.2e ',count,residual,residual/residual_initial,toc(global_timer));
        fprintf('\nIter: %3d \t r: %5.2e \t N(r): %5.2e \t time: %5.2e ',count,residual,residual/residual_initial,toc(global_timer));
        
        %% RESHAPE THE DEFORMATION VECTOR TO TENSOR
        ul_at_n = ul;
        ul(:,1)=ul_list(1:ndim:end,1);
        ul(:,2)=ul_list(2:ndim:end,1);
        if (ndim==3) ul(:,3)=ul_list(3:ndim:end,1);  end
        %% PLOT THE SOLVED SYSTEM UNDER NEWTON-RAPHSON ITERATION
        du= ul-ul_at_n;
        
        if (0);PlotMesh(ul+xl(:,3:3+ndim-1), connect_list,fixed_nodes,disp_nodes); end;
        %% EXIT LOOP IF SATISFIED (NUMBER OF STEPS CAN BE ADDED TOO.)
        % Set the xl vector to xl+ul for the next loading step.
        
        if abs(residual)<tol
            old_dof_list=dof_list;
            old_ul_list=ul_list;
            
            lambda=lambda+dlam;
            x_curr=x_reff+ul;
            x_tec(:,ndim+1:2*ndim)=x_curr;
            fprintf('\n');
            break;
        elseif count==22
            fprintf('\n NO CONVERGENCE, residual= %6.4f \n \n',residual);
            fprintf(file01,'NO CONVERGENCE, residual= %6.4f \n \n',residual);
            stop;
        end
    end
    %% PLOT AT EACH CONVERGED SOLUTION FOR EACH LOAD STEP.
    if (plt_each_load_step); PlotMesh(ul+xl(:,3:3+ndim-1), connect_list,fixed_nodes,disp_nodes); end;
end

toc
%% PLOT FIGURE;
if (plt_final_configuration); PlotMesh(x_curr, connect_list,fixed_nodes,disp_nodes); end;

%% CLOSE THE FILE
fclose(file01);
fclose(file02);
toc(global_timer)
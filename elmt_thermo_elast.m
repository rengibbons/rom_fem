function [K_el,Fd,history_1] = elmt_thermo_elast(sys_vars,ul,old_ul,x_reff,history_0,history_1)
% MAIN PURPOSE: GET STIFFNESS AND RESIDUAL VALUE (s,r)

% _________________ ELMT_THERMO_ELAST.m ____________________

% System parameters
debug     = sys_vars.debug;
file01    = sys_vars.file01;
node_p_el = sys_vars.node_p_el;
n_dim     = sys_vars.n_dim;
dt        = sys_vars.t_para.dt;

% Material parameters
k   = sys_vars.mat_para.k;
c   = sys_vars.mat_para.c;
rho = sys_vars.mat_para.rho;
T0  = sys_vars.mat_para.T0;

ndof_p_node = n_dim + 1;

% Get quadrature information - EVALUATE POSITIONS AND WEIGHTS
if (debug)
    fprintf(file01,'\n GAUSS QUADRATURE - elmt_thermo_elast.m \n ')
end

[~,int_points,int_weights] = gauss_int(sys_vars);
intp_p_node = size(int_points,2);

% Compute current geometry (n_dim+1 x nnode)
x_curr = [x_reff(:,1:n_dim)+ul(:,1:n_dim), ul(:,n_dim+1)];

if (debug)
    mprint(x_reff','REFERENCE COORDINATES (X,Y,Z)',file01);
    mprint(ul','Deformation vector (X,Y,Z)',file01);
    mprint(x_curr','SPATIAL COORDINATES (X,Y,Z)',file01);
end

% Initialize stiffness matrix and residual components
K_el  = zeros(node_p_el*ndof_p_node,node_p_el*ndof_p_node);
K_geo = zeros(node_p_el,node_p_el);

Fd = zeros(node_p_el*ndof_p_node,1);

%% Parse out element dofs
u_dof = zeros(n_dim*node_p_el,1);
T_dof = zeros(node_p_el,1);
for inode = 1 : node_p_el
    if n_dim==2
        u_dof(n_dim*inode-1:n_dim*inode) = [ndof_p_node*inode-2 ndof_p_node*inode-1]';
        T_dof(inode)                   = ndof_p_node*inode;
    elseif n_dim==3
        % TODO: dof assignments for n_dim=3
    end
end

%  ---------------------------------------------------------------------
%
%      LOOP OVER THE INTGRATION POINTS
%
%  ---------------------------------------------------------------------
%%
for intp = 1 : intp_p_node
    %% Compute geometric properties
    % Evaluate shape functions & derivatives
    [N_shape, dNdxi_shape] = shape_fun(int_points(:,intp),sys_vars);
    
    % Compute Jacobian: REFERENCE CONFIGURATION
    [jac_ref,det_jac_ref] = jacobian_mat(x_reff(:,1:n_dim)',dNdxi_shape,sys_vars);
    
    % Compute Jacobian: CURRENT CONFIGURATION
    [jac_cur,det_jac_cur] = jacobian_mat(x_curr(:,1:n_dim)',dNdxi_shape,sys_vars);
    
    % Compute deformation gradient F = je*Je^-1
    F = jac_cur * inv(jac_ref);
    J = det_jac_cur / det_jac_ref;
    
    if n_dim==2
        F33 = zeros(3);
        F33(1:n_dim,1:n_dim) = F;
        F33(n_dim+1,n_dim+1) = 1;
        F = F33;
    end

    % Compute infinitesimal volume element
    dv_cur = det_jac_cur * int_weights(intp);
    
    %% Compute B-operator
    %       B-operator looks like that
    %       | Nk,x    0       0|
    %       |   0    Nk,y     0|
    %       |   0     0    Nk,z|
    %       | Nk,y   Nk,x     0|
    %       |    0   Nk,z  Nk,y|
    %       | Nk,z    0    Nk,x|
    dNidx = inv(jac_cur)'*dNdxi_shape';
   
    if (n_dim==3)
        B_ul = zeros(6,3*node_p_el);
        B_ul(1,1:n_dim:n_dim*node_p_el) = dNidx(1,:);
        B_ul(2,2:n_dim:n_dim*node_p_el) = dNidx(2,:);
        B_ul(3,3:n_dim:n_dim*node_p_el) = dNidx(3,:);
        
        B_ul(4,1:n_dim:n_dim*node_p_el) = dNidx(2,:);
        B_ul(4,2:n_dim:n_dim*node_p_el) = dNidx(1,:);
        
        B_ul(5,2:n_dim:n_dim*node_p_el) = dNidx(3,:);
        B_ul(5,3:n_dim:n_dim*node_p_el) = dNidx(2,:);
        
        B_ul(6,1:n_dim:n_dim*node_p_el) = dNidx(3,:);
        B_ul(6,3:n_dim:n_dim*node_p_el) = dNidx(1,:);
    elseif (n_dim==2)
        B_ul = zeros(3,2*node_p_el);
        B_ul(1,1:n_dim:n_dim*node_p_el) = dNidx(1,:);
        B_ul(2,2:n_dim:n_dim*node_p_el) = dNidx(2,:);
        
        B_ul(3,1:n_dim:n_dim*node_p_el) = dNidx(2,:);
        B_ul(3,2:n_dim:n_dim*node_p_el) = dNidx(1,:);
    end

    %% Compute temperature field and gradients
    T_new  = N_shape *     (ul(:,n_dim+1)+T0);
    T_old  = N_shape * (old_ul(:,n_dim+1)+T0);
    grad_T = dNidx   *     (ul(:,n_dim+1)+T0);

    %% Compute material model
    [cto_mech_voigt,cto_ther_voigt,sigma_voigt,sigma] = mat_thermo_elast_neoh(sys_vars,F,J,T_new,debug,file01);

    %% Compute stiffness matrix components
    % Compute material stiffness matrix
    K_el(u_dof,u_dof) = K_el(u_dof,u_dof) + B_ul'*cto_mech_voigt*B_ul*dv_cur;

    % Compute geometric stiffness matrix
    K_geo = K_geo + dNidx'*sigma*dNidx*dv_cur;

    % Compute coupled mechanical-thermal stiffness matrix
    K_el(u_dof,T_dof) = K_el(u_dof,T_dof) + B_ul'*cto_ther_voigt*N_shape*dv_cur;

    % Compute thermal stiffness matrix
    K_el(T_dof,T_dof) = K_el(T_dof,T_dof) + J*(dNidx'*k*dNidx + N_shape'*rho*c/dt*N_shape)*dv_cur;
    
    %% Compute residuals
    % Compute mechanical residual
    Fd(u_dof) = Fd(u_dof) - B_ul'*sigma_voigt*dv_cur;

    if (sys_vars.heat_src.heat_src_on && ...
       (strcmp(sys_vars.model,'modelI') || strcmp(sys_vars.model,'modelII')) && ...
       history_1.it<ceil(sys_vars.n_time_step/2)) ...%sys_vars.time_step<ceil(sys_vars.n_time_step/2)) ...
       || ...
       (sys_vars.heat_src.heat_src_on && ...
       strcmp(sys_vars.model,'modelIII'))
        % Compute heat source
%         Q = compute_heat_source(sys_vars,N_shape,history_1.x_curr_fixed_el,history_1.it); old
        Q = compute_heat_source(sys_vars,N_shape,history_1);
    else
        Q = 0;
    end
%     Q = compute_heat_source(sys_vars,N_shape,history_1);
    
    % Compute thermal residual
    Fd(T_dof) = Fd(T_dof) - J * k*dNidx'*grad_T                 * dv_cur ...
                          - J * rho*c*N_shape'*(T_new-T_old)/dt * dv_cur ...
                          + J * rho*N_shape'*Q                  * dv_cur;

end % loop over gp

%% Assign geometric stiffness components
for ii = 1 : node_p_el
    for jj = 1 : node_p_el
        ii_ind = (n_dim+1)*ii-n_dim : (n_dim+1)*ii-1;
        jj_ind = (n_dim+1)*jj-n_dim : (n_dim+1)*jj-1;
        K_el(ii_ind,jj_ind) = K_el(ii_ind,jj_ind) + K_geo(ii,jj) * eye(n_dim);
    end
end

end

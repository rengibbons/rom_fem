function [K_el,Fd] = elmt_elast_neoh(sys_vars,ul,x_reff)

% MAIN PURPOSE: GET STIFFNESS AND RESIDUAL VALUE (s,r)

% _________________ ELMT_ELAST_NEOH.m ____________________

% System parameters
debug     = sys_vars.debug;
file01    = sys_vars.file01;
node_p_el = sys_vars.node_p_el;
n_dim     = sys_vars.n_dim;

ndof_p_node = n_dim;

% Get quadrature information - EVALUATE POSITIONS AND WEIGHTS
if (debug)
    fprintf(file01,'\n GAUSS QUADRATURE - elmt_elast_neoh.m \n ')
end

[~,int_points,int_weights] = gauss_int(sys_vars);
intp_p_node = size(int_points,2);

% Compute current geometry (n_dim x nnode)
x_curr = x_reff + ul;

if (debug)
    mprint(x_reff','REFERENCE COORDINATES (X,Y,Z)',file01);
    mprint(ul','Deformation vector (X,Y,Z)',file01);
    mprint(x_curr','SPATIAL COORDINATES (X,Y,Z)',file01);
end

% Initialize stiffness matrix and residual components
K_el  = zeros(node_p_el*ndof_p_node,node_p_el*ndof_p_node);
K_geo = zeros(node_p_el,node_p_el);

Fd = zeros(node_p_el*ndof_p_node,1);
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
    [jac_ref,det_jac_ref] = jacobian_mat(x_reff',dNdxi_shape,sys_vars);
    
    % Compute Jacobian: CURRENT CONFIGURATION
    [jac_cur,det_jac_cur] = jacobian_mat(x_curr',dNdxi_shape,sys_vars);
    
    % Compute deformation gradient F = je*Je^-1
    F = jac_cur * inv(jac_ref);
    J = det_jac_cur / det_jac_ref;
    
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

    %% Compute material model
    [cto_mech_voigt,sigma_voigt,sigma] = mat_elast_neoh(sys_vars,F,J,debug,file01);
    
    %% Compute stiffness matrix components
    % Compute material stiffness matrix
    K_el = K_el + B_ul'*cto_mech_voigt*B_ul*dv_cur;

    % Compute geometric stiffness matrix
    K_geo = K_geo + dNidx'*sigma*dNidx*dv_cur;
    
    %% Compute mechanical residual
    Fd = Fd - B_ul'*sigma_voigt'*dv_cur;

end % loop over gp

%% Assign geometric stiffness components
for ii = 1 : node_p_el
    for jj = 1 : node_p_el
        ii_ind = (ii-1)*n_dim+1 : n_dim*ii;
        jj_ind = (jj-1)*n_dim+1 : n_dim*jj;
        K_el(ii_ind,jj_ind) = K_el(ii_ind,jj_ind) + K_geo(ii,jj) * eye(n_dim);
    end
end

end
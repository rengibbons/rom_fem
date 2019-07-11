function [k_barbar,F_con,add_field] = mixed_elmt13(d,ul,du,xl,add_field)
% MAIN PURPOSE: GET STIFFNESS AND RESIDUAL VALUE (s,r)
% s=zeros(total_node_no,total_node_no); r=zeros(total_node_no,1);
% mat_para=[8,392,0]; % [g12,lambda,gc]
% mat_para=[8,392,0]; % [g12,lambda,gc]

% _________________ MIXED_ELMT13.m ____________________
global node_p_el int_rule debug file01 total_node_no ndim element_type node_p_elm

nodexndim=node_p_elm*ndim;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
du_list=zeros(nodexndim,1);
du_list(1:ndim:end,1)=du(:,1);
du_list(2:ndim:end,1)=du(:,2);
if (ndim==3) du_list(3:ndim:end,1)=du(:,3); end

%%%%%%%%%%%%%%%%%%     I N I T I L I  Z E

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K_mat=zeros(nodexndim,nodexndim);
K_geo=zeros(nodexndim,nodexndim);
K_geo_reduced=zeros(node_p_elm,node_p_elm);
Fd=zeros(nodexndim,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k_uu = zeros(nodexndim,nodexndim);
k_up = zeros(nodexndim,1);
k_Jp = zeros(1,1);
k_JJ = zeros(1,1);

F_p=0;
F_J=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% SHAPE FUNCTIONS FOR J AND p FIELD.
N_J=1;
N_p=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%      COMPUTE THE PRESSURE FIELD OR DILITATION FIELD
J_tilde  = add_field(1,1);
p_tilde  = add_field(2,1);
inv_k_Jp = add_field(3,1); % one-by-one
k_up = add_field(4:nodexndim+3,1);
k_JJ = add_field(nodexndim+4,1);
F_p  = add_field(nodexndim+5,1);
F_j  = add_field(nodexndim+6,1);

d_J = inv_k_Jp'* (F_p - k_up'*du_list);  % I need to get du
d_p = inv_k_Jp'* (F_j - k_JJ*d_J);

%fprintf('\n J = %f , p = %f \n',J_tilde,p_tilde)

p_tilde = p_tilde + d_p;
J_tilde = J_tilde + d_J;


%fprintf('UPDATED: J = %f , p = %f \n',J_tilde,p_tilde)
%fprintf('WITH: d_J = %f , d_p = %f \n',d_J,d_p)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%             R E - I N I T I L I  Z E
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k_up = zeros(nodexndim,1);
k_Jp = zeros(1,1);
k_JJ = zeros(1,1);

F_p=0;
F_J=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     Get quadrature information - EVALUATE POSITIONS AND W8S
if (debug); fprintf(file01,'\n GAUSS QUADRATURE - elmt13.m \n '); end;
[int_points,int_weights]=gauss_int(int_rule,debug,file01);

%     Compute current geometry (ndimXnnode)
x_reff=xl;
x_curr=x_reff+ul;

if (debug)
    mprint(x_reff','REFERENCE COORDINATES (X,Y,Z)',file01);
    mprint(ul','Deformation vector (X,Y,Z)',file01);
    mprint(x_curr','SPATIAL COORDINATES (X,Y,Z)',file01);
end

%  ---------------------------------------------------------------------
%
%      LOOP OVER THE INTGRATION POINTS
%
%  ---------------------------------------------------------------------



for intip=1:size(int_points,2)
    if (debug); fprintf(file01,'\n (%d)^th LOOP OVER THE INTGRATION POINTS - elmt13.m  _____ \n ',intip); end;
    
    % Evaluate shape functions & derivatives
    [N_shape, dNdxi_shape]=shape_fun(int_points(:,intip),element_type,debug,file01);
    
    if (debug); fprintf(file01,'\nREFERENCE'); end;
    % Compute Jacobian: REFERENCE CONFIGURATION
    [jac_ref,det_jac_ref]=jacobian_mat(x_reff',dNdxi_shape,debug,file01);
    
    if (debug); fprintf(file01,'\nCURRENT'); end;
    % Compute Jacobian: CURRENT CONFIGURATION
    [jac_cur,det_jac_cur]=jacobian_mat(x_curr',dNdxi_shape,debug,file01);
    
    % Compute deformation gradient F = je*Je^-1
    F=jac_cur*inv(jac_ref);
    %     F=jac_ref\jac_cur; ! NOT THE SAME THING
    det_F=det_jac_cur/det_jac_ref;
    
    if (debug)
        mprint(jac_ref,'jac_ref ',file01);
        mprint(jac_cur,'jac_cur ',file01);
        mprint(F,'DEFORMATION GRADIENT',file01);
        mprint(det_F,'DETERMINANT OF DEFORMATION GRADIENT',file01);
        mprint(det_jac_ref,'det_jac_ref ',file01);
        mprint(det_jac_cur,'det_jac_cur ',file01);
    end
    
    if(debug)
        mprint(F,'DEFORMATION GRADIENT',file01);
        mprint(det_F,'DETERMINANT OF DEFORMATION GRADIENT',file01);
    end
    
    % Compute B-operator
    %       B-operator looks like that
    %       | Nk,x    0       0|
    %       |   0    Nk,y     0|
    %       |   0     0    Nk,z|
    %       | Nk,y   Nk,x     0|
    %       |    0   Nk,z  Nk,y|
    %       | Nk,z    0    Nk,x|
    
    inv_jac_cur=inv(jac_cur);
    dNidx=inv_jac_cur'*dNdxi_shape';
    
    
    if (ndim==3)
        B=zeros(6,3*node_p_elm);
        B(1,1:ndim:ndim*node_p_elm)=dNidx(1,:);
        B(2,2:ndim:ndim*node_p_elm)=dNidx(2,:);
        B(3,3:ndim:ndim*node_p_elm)=dNidx(3,:);
        
        B(4,1:ndim:ndim*node_p_elm)=dNidx(2,:);
        B(4,2:ndim:ndim*node_p_elm)=dNidx(1,:);
        
        B(5,2:ndim:ndim*node_p_elm)=dNidx(3,:);
        B(5,3:ndim:ndim*node_p_elm)=dNidx(2,:);
        
        B(6,1:ndim:ndim*node_p_elm)=dNidx(3,:);
        B(6,3:ndim:ndim*node_p_elm)=dNidx(1,:);
        
    elseif (ndim==2)
        B=zeros(3,2*node_p_elm);
        
        B(1,1:ndim:ndim*node_p_elm)=dNidx(1,:);
        B(2,2:ndim:ndim*node_p_elm)=dNidx(2,:);
        
        B(3,1:ndim:ndim*node_p_elm)=dNidx(2,:);
        B(3,2:ndim:ndim*node_p_elm)=dNidx(1,:);
        
    end
    
    if (debug)
        mprint(B,'B OPERATOR',file01);
        mprint(dNidx,'DERIVATIVES OF SHAPE FUNCTIONS WITH RESPECT TO CARTESIAN COORDINATES',file01);
    end
    
    % compute infinitesimal volume element
    dv_cur = det_jac_cur*int_weights(intip);
    dv_ref = det_jac_ref*int_weights(intip);
    %%% MATERIAL MODEL
    
    %         [cons_tangent_voigt,cauchy_sigma_voigt,cauchy_sigma]=mate01(F,det_F,d,debug,file01);
    %         [cons_tangent_voigt,cauchy_sigma_voigt,cauchy_sigma]=mate01_2d(F,det_F,d,debug,file01);
    %         [cons_tangent_voigt,cauchy_sigma_voigt,cauchy_sigma]=mate01_growing(F,det_F,d,debug,file01);
%         [cons_tangent_voigt,cauchy_sigma_voigt,cauchy_sigma,J_tilde, p_tilde, dPsi_vol_dJ_tilde, d2Psi_vol_dJ_tilde2]=mate03d_vol_iso(F,det_F,d,J_tilde,p_tilde,debug,file01);
 [cons_tangent_voigt,cauchy_sigma_voigt,cauchy_sigma,J_tilde, p_tilde, dPsi_vol_dJ_tilde, d2Psi_vol_dJ_tilde2,det_Fe]=mate03d_vol_iso_growing(F,det_F,d,J_tilde,p_tilde,debug,file01);
%     [cons_tangent_voigt,cauchy_sigma_voigt,cauchy_sigma,J_tilde, p_tilde, dPsi_vol_dJ_tilde, d2Psi_vol_dJ_tilde2]=mate03d_vol_iso_2d(F,det_F,d,J_tilde,p_tilde,debug,file01);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%% COMPUTE K_uu (24X24)
    
    % Compute MATerial STIFFNESS matrix
    K_mat=K_mat+B'*cons_tangent_voigt*B*dv_cur;
    % Compute GEOmetric STIFFNESS matrix
    K_geo_reduced=K_geo_reduced+dNidx'*cauchy_sigma*dNidx*dv_cur; %'
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%% COMPUTE K_pu (24X1) -> when pressure has the shape function -1-.
    if (ndim==2) m = [1 1 0]'; end
    if (ndim==3) m = [1 1 1 0 0 0]'; end
    
    k_up=k_up+(B'*m).*N_p.*det_Fe.*dv_ref; %'
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%% COMPUTE K_jp (1X1) -> when pressure has the shape function -1-.
    k_Jp=k_Jp-N_J'*N_p*dv_ref;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%% COMPUTE K_jj (1X1) -> when pressure has the shape function -1-.
    
    k_JJ = k_JJ + N_J'* d2Psi_vol_dJ_tilde2 * N_J* dv_ref; %'
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Compute RESIDUAL
    %%%  F_u   24x1
    Fd=Fd-B'*cauchy_sigma_voigt'*dv_cur;
    
    %%%  F_p   1x1         % det_F -> det_Fe
     F_p = F_p - N_p'*(det_Fe - J_tilde)*dv_ref; %'
    
    %%%  F_J   1x1
    F_J = F_J - N_J'*(dPsi_vol_dJ_tilde - p_tilde)*dv_ref; %'
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

if (ndim==3)
    for tt=0:node_p_elm-1
        K_geo(1+3*tt,1:3:3*node_p_elm)=K_geo_reduced(tt+1,:);
        K_geo(2+3*tt,2:3:3*node_p_elm)=K_geo_reduced(tt+1,:);
        K_geo(3+3*tt,3:3:3*node_p_elm)=K_geo_reduced(tt+1,:);
    end
elseif (ndim==2)
    for tt=0:node_p_elm-1
        K_geo(1+2*tt,1:2:2*node_p_elm)=K_geo_reduced(tt+1,:);
        K_geo(2+2*tt,2:2:2*node_p_elm)=K_geo_reduced(tt+1,:);
    end
end

k_uu=K_mat+K_geo;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%  STATIC CONDENSATION

%%%%%%%%%%%%%%%%%%%%%  COMPUTE left HAND SIDE K_barbar d_u = F_con
inv_k_Jp = inv(k_Jp);

% Compute K_bar
k_bar= inv_k_Jp*k_JJ*inv_k_Jp' ; %'
k_barbar= k_uu + k_up*k_bar*k_up'; %'

%%%%%%%%%%%%%%%%%%%%%  COMPUTE RIGHT HAND SIDE

F_con = Fd + k_up*k_bar*F_p-k_up*inv_k_Jp*F_J;


% Residuals
% Res_j = det_F - J_tilde
% Res_p = dPsi_vol_dJ_tilde - p_tilde


% update add_field
%%%%%%      COMPUTE THE PRESSURE FIELD OR DILITATION FIELD
add_field(1,1)= J_tilde;
add_field(2,1)= p_tilde;
add_field(3,1)= inv_k_Jp ;
add_field(4:nodexndim+3,1)=k_up;
add_field(nodexndim+4,1)= k_JJ;
add_field(nodexndim+5,1)= F_p;
add_field(nodexndim+6,1)= F_j;


end
%% INPUT FILE for FES software

mesh_input         = 'patch_test';
reference_geometry = 'quad'; % (2d) quad (3d) cube , circular, complex
element_type       = 'QUAD4';
physics            = 'elastic_neoh';
n_dim              = 2;
node_p_el          = 4;  
int_rule           = '2x2';

n_time_step = 1;
dt          = 1;
tmax        = n_time_step * dt;

E     = 100;
nu    = 0.;%4999;
G     = E / (2 * (1 + nu));
kappa = E / (3 * (1 - 2*nu));

%% EDGE BOUNDARY
% dir,v,i,j,k -> 3d
% dir,v,i,j -> 2d
% % uncomment this for standard footing problem
% define_BC.EBOUN=[2 0.0 1 1];
% define_BC.EBOUN=[2 0.0 0 1;
%                  1 0.0 1 0];

%% COORDINATE BOUNDARY
% x,y,z,i,j,k -> 3D
% x,y,i,j -> 2D
if 1
    define_BC.CBOUN=[0.0 0.0 1 1];
else
    define_BC.CBOUN=[0.0 0.0 1 1;
                     1.0 0.0 1 1];
end
%% COORDINATE FORCE
% x,y,z,i,j,k -> 3D
% x,y,i,j -> 2D
% define_BC.CFORC=[0.00 1.0 0. -2];

%% COORDINATE DISP
% define_BC.CDISP=[0.0 1.0 0.0 0.1;
%                  1.0 1.0 -0.1 0.0];

define_BC.CDISP=[0.4000         0    2.0000   -0.4000;
                 1.0000         0    5.0000   -1.0000;
                      0    0.4000    0.4000         0;
                 1.0000    0.6000    5.6000   -1.0000;
                      0    1.0000    1.0000         0;
                 0.5000    1.0000    3.5000   -0.5000;
                 1.0000    1.0000    6.0000   -1.0000];
             
%% EDGE DISP
% dir,v,dx,dy,dz -> 3d
% dir,v,dx,dy -> 2d
% define_BC.EDISP=[2 1 0 0.1];

t_para = struct('dt',dt,...
                'tmax',tmax);

mat_para = struct('G',G,...
                  'kappa',kappa);
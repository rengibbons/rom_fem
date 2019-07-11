%% INPUT FILE for FES software

mesh_input         = 'square04x04';
reference_geometry = 'quad'; % (2d) quad (3d) cube , circular, complex
element_type       = 'QUAD4';
physics            = 'elastic_neoh';
n_dim              = 2;
node_p_el          = 4;  
int_rule           = '2x2';

disp  = 0.02;

n_time_step = 10;
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
define_BC.EBOUN=[2 0.0 1 1];
% define_BC.EBOUN=[2 0.0 0 1;
%                  1 0.0 1 0];

%% COORDINATE BOUNDARY
% x,y,z,i,j,k -> 3D
% x,y,i,j -> 2D
% define_BC.CBOUN=[0.0 0.0 1 1; 
%                  1.0 0.0 1 1];

%% COORDINATE FORCE
% x,y,z,i,j,k -> 3D
% x,y,i,j -> 2D
% define_BC.CFORC=[0.00 1.0 0. -2];

%% COORDINATE DISP
% define_BC.CDISP=[0.00 1.0 0.0 disp];

%% EDGE DISP
% dir,v,dx,dy,dz -> 3d
% dir,v,dx,dy -> 2d
define_BC.EDISP=[2 0.1 0 disp];

t_para = struct('dt',dt,...
                'tmax',tmax);

mat_para = struct('E',E,...
                  'nu',nu,...
                  'G',G,...
                  'kappa',kappa);

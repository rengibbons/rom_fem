%% INPUT FILE for FES software

mesh_input         = 'footing04x04';
% mesh_input         = 'footing20x20';
reference_geometry = 'quad'; % (2d) quad (3d) cube , circular, complex
element_type       = 'QUAD4';
physics            = 'elastic';
ndim               = 2;
node_p_el          = 4;  
int_rule           = '2x2';

width = 0.4;
disp  = -0.13;

tmax = 5;
dt   = 1;

E      = 100;
nu     = 0.3;
G      = E / 2 / (1+nu);
lambda = E * nu / (1+nu) / (1-2*nu);

mat_para = struct('G',nu,...
                  'lambda',lambda);

%% EDGE BOUNDARY
% dir,v,i,j,k -> 3d
% dir,v,i,j -> 2d
% % uncomment this for standard footing problem
define_BC.EBOUN=[2 0.0 1 1];
%                  1 0.0 1 0;
%                  1 1.0 1 1];

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
define_BC.EDISP=[2 1.0 0 disp];

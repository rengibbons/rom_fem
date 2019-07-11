%% INPUT FILE for FES software

mesh_input         = 'footing04x04';
% mesh_input         = 'footing20x20';
reference_geometry = 'quad'; % (2d) quad (3d) cube , circular, complex
element_type       = 'QUAD4';
physics            = 'thermo_elastic';
ndim               = 2;
node_p_el          = 4;  
int_rule           = '2x2';

width = 0.4;
disp  = -0.13;

tmax = 5;
dt   = 1;

% Mechanical properties
E      = 100;
nu     = 0.3;
kappa = E / 3 / (1 - 2*nu);
G     = E / 2 / (1 + nu);

% Thermal properties
alpha = 17.3e-5; % thermal expansion
k     = 16.2e-6; % conductivity
c     = 3.756;   % capacity
T0    = 20;      % initial temp

mat_para = struct('kappa',kappa,...
                  'G',G,...
                  'T0',T0,...
                  'alpha',alpha,...
                  'k',k,...
                  'c',c);

%% EDGE BOUNDARY
% dir,v,i,j,k -> 3d
% dir,v,i,j   -> 2d
% % dir,v,i,j,k,T -> 3d (thermoelastic)
% % dir,v,i,j,T   -> 2d (thermoelastic)
% % uncomment this for standard footing problem
define_BC.EBOUN=[2 0.0 1 1 1];
%                  1 0.0 1 0 0;
%                  1 1.0 1 1 0];       

%% COORDINATE BOUNDARY
% x,y,z,i,j,k -> 3D
% x,y,i,j     -> 2D
% x,y,z,i,j,k,T -> 3D (thermoelastic) 
% x,y,i,j,T     -> 2D (thermoelastic)
% define_BC.CBOUN=[0.0 0.0 1 1 1;
%                  0.5 0.0 0 0 1;
%                  1.0 0.0 1 1 1];

%% COORDINATE FORCE
% x,y,z,i,j,k -> 3D
% x,y,i,j     -> 2D
% x,y,z,i,j,k,flux -> 3D (thermoelastic)
% x,y,i,j,flux     -> 2D (thermoelastic)
% define_BC.CFORC=[0.00 1.0 0. -2];

%% EDGE DISP
% dir,v,dx,dy,dz -> 3d
% dir,v,dx,dy    -> 2d
% dir,v,dx,dy,dz,dT -> 3d (thermoelastic)
% dir,v,dx,dy,dT    -> 2d (thermoelastic)
% define_BC.EDISP=[2 1.0 0 disp 17];

%% COORDINATE DISP
% x,y,z,dx,dy,dz -> 3d
% x,y,dx,dy      -> 2d
% x,y,z,dx,dy,dz,dT -> 3d (thermoelastic)
% x,y,dx,dy,dT      -> 2d (thermoelastic)
% define_BC.CDISP=[0.00 1.0 0.0 disp];

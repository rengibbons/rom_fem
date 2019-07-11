%% INPUT FILE for FES software

mesh_input         = 'square10x10';
reference_geometry = 'quad';
element_type       = 'QUAD4';
physics            = 'thermo_elastic';
n_dim              = 2;
node_p_el          = 4;  
int_rule           = '2x2';

n_time_step = 10;
dt   = 0.1;
tmax = dt*n_time_step;

% Mechanical properties
E     = 100;
nu    = 0.3;
kappa = E / (3 * (1 - 2*nu));
G     = E / (2 * (1 + nu));

% Thermal properties
alpha = 1e-4;%17.3e-3; % thermal expansion
k     = 60.5;    % conductivity
c     = 452;     % capacity
rho   = 7850;    % desnsity
T0    = 20;      % initial temp

% Thermal properties
% alpha = 5e-3;    % thermal expansion
% k     = 1e-8;    % conductivity
% c     = 1e-8;    % capacity
% rho   = 1.1;     % density
% T0    = 20;      % initial temp
              
% Heat source properties
heat_src_on = 1;
s0          = 2e4;
r           = 0.01;

% Beginning and final points of heat source
x0 = 0.9*[0.1,0.1];
xf = 0.1*[0.1,0.1];

x_path = linspace(x0(1),xf(1),n_time_step);
y_path = linspace(x0(2),xf(2),n_time_step);
heat_path = [x_path' y_path'];

%% EDGE BOUNDARY
% dir,v,i,j,k -> 3d
% dir,v,i,j   -> 2d
% dir,v,i,j,k,T -> 3d (thermoelastic)
% dir,v,i,j,T   -> 2d (thermoelastic)
define_BC.EBOUN=[2 0.0 0 1 0;
                 1 0.0 1 0 0];
% define_BC.EBOUN=[2 0.0 0 1 0];


%% EDGE DISP
% dir,v,dx,dy,dz -> 3d
% dir,v,dx,dy    -> 2d
% dir,v,dx,dy,dz,dT -> 3d (thermoelastic)
% dir,v,dx,dy,dT    -> 2d (thermoelastic)
% define_BC.EDISP=[2 1.0 0 disp Tinc];
% define_BC.EDISP=[2 1.0 0 disp 0];
% Tinc = 30;
% define_BC.EDISP=[2 0.1 0 0 Tinc];

%% COORDINATE BOUNDARY
% x,y,z,i,j,k -> 3D
% x,y,i,j     -> 2D
% x,y,z,i,j,k,T -> 3D (thermoelastic) 
% x,y,i,j,T     -> 2D (thermoelastic)
% define_BC.CBOUN=[0.0 0.0 1 1 0];
%                  0.1 0.0 0 1 0];
%                  0.5 0.0 0 0 1;
%                  1.0 0.0 1 1 1];

%% COORDINATE FORCE
% x,y,z,i,j,k -> 3D
% x,y,i,j     -> 2D
% x,y,z,i,j,k,flux -> 3D (thermoelastic)
% x,y,i,j,flux     -> 2D (thermoelastic)
% define_BC.CFORC=[0.00 1.0 0. -2];

%% COORDINATE DISP
% x,y,z,dx,dy,dz -> 3d
% x,y,dx,dy      -> 2d
% x,y,z,dx,dy,dz,dT -> 3d (thermoelastic)
% x,y,dx,dy,dT      -> 2d (thermoelastic)
% define_BC.CDISP=[0.00 1.0 0.0 disp];

t_para = struct('dt',dt,...
                'tmax',tmax);

mat_para = struct('E',E,...
                  'nu',nu,...
                  'kappa',kappa,...
                  'G',G,...
                  'T0',T0,...
                  'alpha',alpha,...
                  'k',k,...
                  'c',c,...
                  'rho',rho);

heat_src = struct('heat_src_on',heat_src_on,...
                  's0',s0,...
                  'r',r,...
                  'heat_path',heat_path);

%% INPUT FILE for FES software

mixed=1;

mesh_input='xx';
reference_geometry='quad'; % (2d) quad (3d) cube , circular, complex
element_type='QUAD4';
ndim=2;
node_p_el=4;  
% int_rule='1point';
int_rule='2x2';


d=zeros(15,1);
d(1:3,1)=[8;392;0.0]; % shear modulus, lambda, growth parameter

% dir,v,i,j,k -> 3d
% dir,v,i,j -> 2d
define_BC.EBOUN=[2 -0.0 1 1];

% x,y,z,i,j,k -> 3D
% x,y,i,j -> 2D
% define_BC.CBOUN=[-.5 -.5 1 1];

% dir,v,dx,dy,dz -> 3d
% dir,v,dx,dy -> 2d
define_BC.EDISP=[2 1 0 1];


        
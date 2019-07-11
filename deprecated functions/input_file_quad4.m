%% INPUT FILE for FES software

% mesh_input='xxf';
% reference_geometry='quad'; % (2d) quad (3d) cube , circular, complex
% element_type='QUAD4';
% ndim=2;
% node_p_el=4;  
% int_rule='2x2';
% 

mesh_input='aa';
reference_geometry='quad'; % (2d) quad (3d) cube , circular, complex
element_type='QUAD4';
ndim=2;
node_p_el=4;  
int_rule='2x2';


d=zeros(15,1);
d(1:3,1)=[304.054;7297.3;0.0]; % shear modulus, lambda, growth parameter

% dir,v,i,j,k -> 3d
% dir,v,i,j -> 2d
%define_BC.EBOUN=[2 -0.5 1 1];
           
             
% x,y,z,i,j,k -> 3D
% x,y,i,j -> 2D
define_BC.CBOUN=[0.0  0.0 1 1; 10. 0. 1 1];
% x,y,z,i,j,k -> 3D
% x,y,i,j -> 2D    1.500625e+01 	  8.749220e-01
define_BC.CFORC=[6.509080e+00 	  1.454864e+00 0. -20.];
%

% dir,v,dx,dy,dz -> 3d
% dir,v,dx,dy -> 2d
% define_BC.EDISP=[2 0.5 0 1];        
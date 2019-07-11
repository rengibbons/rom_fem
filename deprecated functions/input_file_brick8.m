%% INPUT FILE for FES software

mesh_input='bloc4';
reference_geometry='cube'; % cube , circular, complex
element_type='BRICK8';  % BRICK27
ndim=3;
node_p_el=8;  %BRCIK8 ELEMENT
int_rule='2x2x2';


d=zeros(15,1);
d(1:3,1)=[8;392;0.050]; % shear modulus, lambda, growth parameter

% dir,v,i,j -> 2d
define_BC.EBOUN=[3 0 1 1 1
                 3 1 1 1 1];
           
             
% x,y,z,i,j,k -> 3D
% x,y,i,j -> 2D
% define_BC.CBOUN=[-.5 -.5 1 1];

% dir,v,dx,dy,dz -> 3d
% dir,v,dx,dy -> 2d
define_BC.EDISP=[3 1 0 0. 0.];


       
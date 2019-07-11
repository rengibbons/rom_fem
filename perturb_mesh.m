clc, close all
folder = 'input';
mesh_input_file = 'square10x10';
mesh_output_file = 'square_perturbed10x10';


nodes = sprintf('input/%s.nodes',mesh_input_file)
coor_list = dlmread(nodes);

xl = 0;
xr = 0.1;
yl = 0;
yu = 0.1;

nn = size(coor_list,1);
x_pert = zeros(nn,1);
y_pert = zeros(nn,1);
for ii = 1 : nn
    if coor_list(ii,3)~=xl && coor_list(ii,3)~=xr && coor_list(ii,4)~=yl && coor_list(ii,4)~=yu
        coor_list(ii,3) = coor_list(ii,3) + 0.005*(rand - 0.5);
        coor_list(ii,4) = coor_list(ii,4) + 0.005*(rand - 0.5);
    end
end

dlmwrite(strcat(folder,'/',mesh_output_file,'.nodes'),coor_list,'delimiter','\t','precision',10)


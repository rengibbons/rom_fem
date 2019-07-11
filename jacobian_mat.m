function [ jac_mat,det_jac_mat] = jacobian_mat(coords,dNdxi_shape,sys_vars)
%JACOBIAN_MAT Summary of this function goes here
%   Detailed explanation goes here
%   Ref: Hughes pg 119

% coords = [x1 x2 x3 x4
%           y1 y2 y3 y4]

debug  = sys_vars.debug;
file01 = sys_vars.file01;

jac_mat     = coords * dNdxi_shape;
det_jac_mat = det(jac_mat);

if (det_jac_mat<0)
    fprintf('JACOBIAN NEGATIVE\n\n\n\nJACOBIAN NEGATIVE\n');
    fprintf(file01,'JACOBIAN NEGATIVE\n\n\n\nJACOBIAN NEGATIVE\n');
end

if(debug)
    mprint(jac_mat,'JACOBIAN MATRIX',file01);
    mprint(det_jac_mat,'DETERMINANT OF JACOBIAN MATRIX',file01);
end
end


function [ residual_value ] = compute_residual( r,bc )
%COMPUTE_RESIDUAL Summary of this function goes here
%   Detailed explanation goes here
r(bc)=[];
residual_value=norm(r,2);

end


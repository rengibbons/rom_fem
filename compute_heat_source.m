function Q = compute_heat_source(sys_vars,N_shape,history_1)

n_dim    = sys_vars.n_dim;
s0       = sys_vars.heat_src.s0;
r0       = sys_vars.heat_src.r;
heat_loc = sys_vars.heat_src.heat_path(history_1.it,:);
% myheatloc = heat_loc

x_gp = N_shape * history_1.x_curr_fixed_el(:,1:2);

if n_dim == 2
    Q = s0 * exp(-((x_gp(1)-heat_loc(1))^2 + (x_gp(2)-heat_loc(2))^2) / r0^2);
else 
    fprintf('Implement 3D heat source term later\n')
end

end

% function Q = compute_heat_source(sys_vars,N_shape,x_curr,history_1)
% 
% n_dim    = sys_vars.n_dim;
% s0       = sys_vars.heat_src.s0;
% r0       = sys_vars.heat_src.r;
% heat_loc = sys_vars.heat_src.heat_path(history_1,:);
% 
% x_gp = N_shape * x_curr(:,1:2);
% 
% if n_dim == 2
%     Q = s0 * exp(-((x_gp(1)-heat_loc(1))^2 + (x_gp(2)-heat_loc(2))^2) / r0^2);
% else 
%     fprintf('Implement 3D heat source term later\n')
% end
% 
% end


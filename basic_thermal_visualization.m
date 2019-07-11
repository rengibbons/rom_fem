function basic_thermal_visualization(connect_list,x_curr,thermal_solution)

figure()
trisurf(connect_list,x_curr(:,1),x_curr(:,2),thermal_solution,'facecolor','interp')

buffer = 0.01;

thermal_min = min(thermal_solution);
thermal_max = max(thermal_solution);

if thermal_min > 0
    thermal_min = thermal_min * (1-buffer);
else
    thermal_min = thermal_min * (1+buffer);
end

if thermal_max > 0
    thermal_max = thermal_max * (1+buffer);
else
    thermal_max = thermal_max * (1-buffer);
end
    
axis([-1e-6,max(x_curr(:,1)),-1e-6,max(x_curr(:,2)),thermal_min,thermal_max])
view(0,90);
colorbar
caxis([thermal_min,thermal_max])
% print('basic_ther_exp_example','-depsc','-tiff')

end
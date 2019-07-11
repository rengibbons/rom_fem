function plot_loading_trend(t_0,tmax,total_load_step,num_load_step);
% PLT_LOADING_TREND Summary of this function goes here
%   Detailed explanation goes here

A=1; % LOAD AMPLITUDE 
dt=tmax/total_load_step;

x_axis=t_0(1):dt:dt*num_load_step;
y_axis=A.*x_axis./(dt*num_load_step);

x_axis_rest=x_axis(end):dt:tmax;
y_axis_rest=x_axis_rest.*0+A;

x_plt=[x_axis x_axis_rest];
y_plt=[y_axis y_axis_rest];


stem(x_plt,y_plt,'fill','LineWidth',2);
hold on;
line(x_plt,y_plt,'LineWidth',4);
xlim([0 tmax*1.05]);
ylim([0 1.05]);
xlabel('t, time');
ylabel('A, Load factor')

end


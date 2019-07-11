clc, clear all, close all

int_n       = 4;
p1          = sqrt(1/3);
xi  = p1*[-1 -1;
           1 -1;
           1  1;
          -1  1];
       
int_weights = [1 1 1 1];

N = zeros(4,1);
N(1) = (1+xi(1,1))*(1+xi(1,2))/4;
N(2) = (1+xi(2,1))*(1+xi(2,2))/4;
N(3) = (1+xi(3,1))*(1+xi(3,2))/4;
N(4) = (1+xi(4,1))*(1+xi(4,2))/4;

x_curr = [0.0 0.0;
          0.5 0.0;
          0.5 0.5;
          0.0 0.5];

% x_curr = [1.0 0.0;
%           1.5 0.3;
%           1.8 0.8;
%           0.4 0.5];

% gp_ana = 0.5/2 * [(1-sqrt(1/3)) (1-sqrt(1/3));
%                   (1+sqrt(1/3)) (1-sqrt(1/3));
%                   (1+sqrt(1/3)) (1+sqrt(1/3));
%                   (1-sqrt(1/3)) (1+sqrt(1/3))]

x_gp = zeros(4,2);
% y_gp = zeros(4,1);

for ii = 1 : int_n
    if ii==1
        Ni = [N(1) N(2) N(3) N(4)];
    elseif ii==2
        Ni = [N(2) N(3) N(4) N(1)];
    elseif ii==3
        Ni = [N(3) N(4) N(1) N(2)];
    elseif ii==4
        Ni = [N(4) N(1) N(2) N(3)];
    end
    Ni
    x_curr
    x_gp(ii,:) = Ni * x_curr;%(:,1);
%     y_gp(ii) = Ni * x_curr(:,2);
    fprintf('ii = %.0f  x = %.5f  y = %.5f\n',[ii,x_gp(ii,1),x_gp(ii,2)])
end

figure, hold on
plot([x_curr(1,1) x_curr(2,1)] , [x_curr(1,2) x_curr(2,2)], 'b')
plot([x_curr(2,1) x_curr(3,1)] , [x_curr(2,2) x_curr(3,2)], 'b')
plot([x_curr(3,1) x_curr(4,1)] , [x_curr(3,2) x_curr(4,2)], 'b')
plot([x_curr(4,1) x_curr(1,1)] , [x_curr(4,2) x_curr(1,2)], 'b')

scatter(x_gp(:,1),x_gp(:,2),'rx')

axis([min(x_curr(:,1))*(1-0.05) max(x_curr(:,1))*(1+0.05) min(x_curr(:,2))*(1-0.05) max(x_curr(:,2))*(1+0.05)])
              
              
              
              
              
              
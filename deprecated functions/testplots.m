clc, clear all, close all
disp('hello')

% for i = 1 : 5
% %     hold off;figure(1);hold on;view(2);set(0,'DefaultFigureVisible','off')
% %     hold off;figure(1);hold on;%view(2);%set(0,'DefaultFigureVisible','off')
%     x = 1 : 10;
%     y = rand(1,10);
% 
%     figure(1), hold off 
%     plot(x,y)
%     
% 
% end
% hold off


x = [0 0.5 1 0.15 0.5 0.85 -0.15 0.5 1.15]';
y = [0 0 0 0.5 0.6 0.5 1 1.2 1]';
T = [5 5 5 10 10 10 20 20 20]';

n = 1000;
[xi, yi] = meshgrid(...
    linspace(min(x),max(x),n),...
    linspace(min(y),max(y),n));

Ti = griddata(x,y,T, xi,yi);

figure
hold on
surface(xi,yi,Ti,'EdgeColor','none','LineStyle','none','FaceLighting','phong')
scatter(x,y)
hold off



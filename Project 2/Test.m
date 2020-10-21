clear all;
close all;

t = 1;

y = 1/3*((t-2)^3) + 1/2*((5-t)*(t-2)^2) + ((4-t)*(t-2)) + 13/3 - (5-3*t)/2

t = 2;

y = 1/3*((t-2)^3) + 1/2*((5-t)*(t-2)^2) + ((4-t)*(t-2)) + 13/3 - (5-3*t)/2



% n1 = linspace(-1,1,3);
% 
% x = heaviside(n1+1).*(1+n1) - heaviside(n1).*((1+n1) - (1-n1));
% 
% n2 = linspace(1,4,4);
% 
% h = 2*heaviside(n2-1) - 2*heaviside(n2-2) + heaviside(n2-2).*(4-n2);
% 
% y1 = conv(x,h);
% 
% plot(y1, 'g', 'LineWidth', 2);
% hold on;
% plot(x, 'r', 'LineWidth', 2);
% plot(h, 'b', 'LineWidth', 2);
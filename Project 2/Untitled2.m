
clear all;
close all;

a = 0.5;
N = 1000;

n = linspace(0,N, 4*N); % I like the smoother curves

load('lineup.mat');

% This is part a
he = [1 zeros(1,998) a]; % There should only be an impulse at n = 0, 1000


% This is part c

a1 = [1 zeros(1,998) a]; % System coefficients for z
b1 = 1;

d = [1 zeros(1,4000)];

her = filter(1, a1, d);

z = filter(b1, a1, y);
plot(y, 'b', 'DisplayName', 'y');

hold on;

plot(z, 'r', 'DisplayName', 'z');

legend('boxon');

sound(y, 8192);

pause(3);

sound(z, 8192);

hoa = conv(he, her);

figure(2);
plot(he, '--g');
hold on;
plot(her, '.r');
plot(hoa, '-b', 'LineWidth', 2);
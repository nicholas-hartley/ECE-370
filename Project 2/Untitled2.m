clear all;
close all;

a = 0.5;
N = 1000;

n = linspace(0,N, 2*N); % I like the smoother curves

load('lineup.mat');

he = [1 zeros(1,998) a]; % Model dirac delta using comparison

d = [1 zeros(1,4000)];

a1 = [1 a];

b1 = 1;

her = filter(b1, a1, d);

z = filter(b1, a1, y);
plot(y);

hold on;

plot(z, 'r');

% sound(y, 8192);
% 
% pause(3);
% 
% sound(z, 8192);

hoa = conv(her, he);

figure(2);
plot(he, '--g');
hold on;
plot(her, '.r');
plot(hoa);
clear all;
% close all;

T = .005;

n = [0:T:20];

a1 = (T * (1 + T^2)^2 * (1i - T))/((1-1i*T)^3 * (T - 2i) + (1+1i*T)^3 * (T + 2i));
a2 = (T * (1 + T^2)^2 * (1i + T))/((1-1i*T)^3 * (T - 2i) + (1+1i*T)^3 * (T + 2i));

yn = a1*((1+1i*T).^n) + a2*((1-1i*T).^n) - (a1 + a2);

yn = yn .* heaviside(n);

yt = 1 - cos(n);

yt = yt .* heaviside(n);

% figure(1);
% 
% plot(n, yt);
% 
% hold on;
% 
% figure(2);
% 
% plot(n, abs(yn));

a = [-2, (1+T^2)];
b  = [T^2, 0, 0];


y0 = [2, 1];
x0 = [0, 0];

x = heaviside(n);

y = recur(a,b,n,x,x0,y0);

tiledlayout(2,2)
stem(nexttile, n, x)

stem(nexttile, n, y)

xlabel('n')
ylabel('y[n]')

stem(nexttile, n, abs(yn))

stem(nexttile, n, yt)

clear all;
close all;

a = [1.5, 0.5];
b  = [1 0 -1];

y0 = [-2, -1];
x0 = [-1, 4];

n = [0:1:20];

x = ((0.5).^(n-1)).*heaviside(n-1);

y = recur(a,b,n,x,x0,y0);

%tiledlayout(1,2)
% stem(nexttile, n, x)

%stem(nexttile, n, y)
stem(n, y)
xlabel('n')
ylabel('y[n]')
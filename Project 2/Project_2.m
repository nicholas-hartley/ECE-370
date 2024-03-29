%% Problem 1

% Problem 1-a&b
clear all;
close all;

a1 = 1;
b1 = [0.5 1 2];

n = [1:1:4];

x = n.*heaviside(n);

ya = filter(b1,a1,x)

a2 = [1 -0.8];
b2  = [0 1];

n = [1:1:20];

x = heaviside(n);

yb = filter(b2,a2,x)

% Problem 1-c
clear all;
close all;

a = -0.8;
b  = [0 1];

y0 = 0;
x0 = 0;

n = [1:1:20];

x = heaviside(n);

y = recur(a,b,n,x,x0,y0);

stem(n,y)
xlabel('n')
ylabel('y[n]')

% Problem 1-d
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

%% Problem 2

clear all;
close all;

tiledlayout(2,2);

for c = 1:4
    T = 0.5^(c*c);

    n = [0:T:20];

    yt = 1 - cos(n);

    yt = yt .* heaviside(n);

    a = [-2, (1+T^2)];
    b  = [T^2, 0, 0];

    y0 = [0, 0];
    x0 = [0, 0];

    x = heaviside(n);

    yn = recur(a,b,n,x,x0,y0);
    
    nexttile
    
    plot(n, yt, 'b', 'DisplayName', 'y(t)' )
    
    hold on;
    
    plot(n, yn, 'r', 'DisplayName', 'y[n]')
    
    title("T is " + T);
    xlabel('n');
    ylabel('y(t) & y[n]');
    
end

%% Problem 3

clear all;
close all;

nx = 0;
Nx = 5;
Samples = 100;

n = linspace(nx,Nx,Samples);

t = linspace(nx,2*Nx,199); % Would be 2*Nx - (Nx-nx)/samples if not inclusive

x = heaviside(n);

h = n.*heaviside(n);

y = conv(x, h);

plot(t, y)

title('conv(x[n], h[n])');
xlabel('n');
ylabel('y[n]');

%% Problem 4

clear all;
close all;

a = 0.5;
N = 1000;

n = linspace(0,N, N); % I like the smoother curves

load('lineup.mat');

% This is part a
he = [1 zeros(1,999) a]; % There should only be an impulse at n = 0, 1000
% Side note matlab is 1 indexed, so n = 1000 correlates to 1001

plot(he)
% This is part c

a1 = [1 zeros(1,999) a]; % System coefficients for z
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
% plot(her, '.-r', 'LineWidth', 2);
% hold on;
plot(hoa, '-b', 'LineWidth', 2);
% plot(he, '--g', 'LineWidth', 2);
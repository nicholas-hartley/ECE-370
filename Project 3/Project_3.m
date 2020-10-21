%% Problem 1

clear all;
close all;

t = [0:0.001:1];
n = mod(11808942, 41);

A = (rand(1,n).*3)+(rand(1,n).*3i);
omega = (rand(1,n).*pi());

xs = SUMCS(t,A,omega);

tiledlayout(2,2);

nexttile

plot(t, real(xs), 'r')
    
title("Real part of xs v.s. t");
xlabel('t');
ylabel('real(xs)');

nexttile

plot(t, imag(xs), 'r')
    
title("Imaginary part of xs v.s. t");
xlabel('t');
ylabel('imag(xs)');

nexttile

plot(t, abs(xs), 'r')
    
title("Magnitude of xs v.s. t");
xlabel('t');
ylabel('abs(xs)');

nexttile

plot(t, angle(xs), 'r')
    
title("Phase angle of xs v.s. t");
xlabel('t');
ylabel('angle(xs)');



function [xs] = SUMCS(t,A,omega)
    xs = 0;
    for i = 1:length(A)
        % Don't do elementwise multiplication of t so that it can generate
        % different sized vectors given different t inputs
        xs = xs + A(i)*exp(j*omega(i)/2*t);
    end
end

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
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

% It'll be at the bottom so that matlab stops yelling

% function [xs] = SUMCS(t,A,omega)
%     xs = 0;
%     for i = 1:length(A)
%         % Don't do elementwise multiplication of t so that it can generate
%         % different sized vectors given different t inputs
%         xs = xs + A(i)*exp(j*omega(i)/2*t);
%     end
% end

%% Problem 2

%{

All of problem 2 was hand written math to find the Fourier Series expansion
coefficients

%}


%% Problem 3

clear all;
close all;

t = [-5:0.001:5];
D11 = mod(11808942, 11);
D4 = mod(11808942, 4);

T = 2;
W = 1;
K = 23 + D11;

xt = FSWave(t, K, T, W);

%% Single xt

% tiledlayout(2,1);
% 
% nexttile
% 
% plot(t, real(xt), '.-r', 'LineWidth', 1.5);
% 
% title("Real part of xt v.s. t");
% xlabel('t');
% ylabel('real(xt)');
% 
% nexttile
% 
% plot(t, imag(xt), '.-b', 'LineWidth', 1.5)
%     
% title("Imaginary part of xt v.s. t");
% xlabel('t');
% ylabel('imag(xt)');

%% Five Plots

tiledlayout(3,2);

nexttile

K = 1 + D4;
xt = FSWave(t, K, T, W);

plot(t, real(xt),  '.-r', 'LineWidth', 1.25)
    
title("Real part of xt v.s. t for K=" + K);
xlabel('t');
ylabel('real(xt)');

nexttile

K = 7 + D4;
xt = FSWave(t, K, T, W);

plot(t, real(xt), '.-r', 'LineWidth', 1.25)
    
title("Real part of xt v.s. t for K=" + K);
xlabel('t');
ylabel('real(xt)');

nexttile

K = 16 + D4;
xt = FSWave(t, K, T, W);

plot(t, real(xt),  '.-r', 'LineWidth', 1.25)
    
title("Real part of xt v.s. t for K=" + K);
xlabel('t');
ylabel('real(xt)');

nexttile

K = 100 + D4;
xt = FSWave(t, K, T, W);

plot(t, real(xt),  '.-r', 'LineWidth', 1.25)
    
title("Real part of xt v.s. t for K=" + K);
xlabel('t');
ylabel('real(xt)');

nexttile

K = 200 + D4;
xt = FSWave(t, K, T, W);

plot(t, real(xt),  '.-r', 'LineWidth', 1.25)
    
title("Real part of xt v.s. t for K=" + K);
xlabel('t');
ylabel('real(xt)');

% Commented out so matlab will stop yelling at me

% function [xs] = SUMCS(t,A,omega)
%     xs = 0;
%     for i = 1:length(A)
%         xs = xs + A(i)*exp(j*omega(i)/2*t);
%     end
% end

% function [xt] = FSWave(t,K,T,W)
%     syms x k
%     expr = (1-4*x^2)*exp(-1i*2*pi*k*x/T);
%     xt = zeros(1,length(t));
%     for i = -K:K % Not assuming discrete time
%         expr2 = subs(expr, k, i);
%         Xk = (1/T)*vpaintegral(expr2, x, -W/4, W/4);
%         % Summing vectors of same size
%         xt = xt + SUMCS(t, Xk, 4*pi()/T*i);
%     end
% end

%% Problem 4

function [xs] = SUMCS(t,A,omega)
    xs = 0;
    for i = 1:length(A)
        xs = xs + A(i)*exp(j*omega(i)/2*t);
    end
end

function [xt] = FSWave(t,K,T,W)
    syms x k
    expr = (1-4*x^2)*exp(-1i*2*pi*k*x/T);
    xt = zeros(1,length(t));
    for i = -K:K % Not assuming discrete time
        expr2 = subs(expr, k, i);
        Xk = (1/T)*vpaintegral(expr2, x, -W/4, W/4);
        % Summing vectors of same size
        xt = xt + SUMCS(t, Xk, 4*pi()/T*i);
    end
end
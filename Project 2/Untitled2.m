clear all;
close all;

t = [-5:0.001:5];
D11 = mod(11808942, 11);
D4 = mod(11808942, 4);

T = 2;
W = 1;
K = 23 + D11;

%% Single xt

% xt = FSWave(t,K,T,W);
% 
% tiledlayout(2,1);
% 
% nexttile
% 
% plot(t, real(xt), 'r')
%     
% title("Real part of xt v.s. t");
% xlabel('t');
% ylabel('real(xt)');
% 
% nexttile
% 
% plot(t, imag(xt), 'r')
%     
% title("Imaginary part of xt v.s. t");
% xlabel('t');
% ylabel('imag(xt)');

%% Five Plots

tiledlayout(3,2);

nexttile

K = 1 + D4;
xt = FSWave(t, K, T, W);

plot(t, real(xt), 'r')
    
title("Real part of xt v.s. t for K=" + K);
xlabel('t');
ylabel('real(xt)');

nexttile

K = 7 + D4;
xt = FSWave(t, K, T, W);

plot(t, real(xt), 'r')
    
title("Real part of xt v.s. t for K=" + K);
xlabel('t');
ylabel('real(xt)');

nexttile

K = 16 + D4;
xt = FSWave(t, K, T, W);

plot(t, real(xt), 'r')
    
title("Real part of xt v.s. t for K=" + K);
xlabel('t');
ylabel('real(xt)');

nexttile

K = 100 + D4;
xt = FSWave(t, K, T, W);

plot(t, real(xt), 'r')
    
title("Real part of xt v.s. t for K=" + K);
xlabel('t');
ylabel('real(xt)');

nexttile

K = 200 + D4;
xt = FSWave(t, K, T, W);

plot(t, real(xt), 'r')
    
title("Real part of xt v.s. t for K=" + K);
xlabel('t');
ylabel('real(xt)');

function [xs] = SUMCS(t,A,omega)
    xs = 0;
    for i = 1:length(A)
        xs = xs + A(i)*exp(j*omega(i)/2*t);
    end
end

function [xt] = FSWave(t,K,T,W)
    syms x
    expr = 1-4*x^2;
    xt = 0;
    for i = -K:K % Not assuming discrete time
        Xk = (1/T)*vpaintegral(expr, x, -W/4, W/4);
        % Summing vectors of same size
        xt = xt + SUMCS(t, Xk, 4*pi()/T*i);
    end
end
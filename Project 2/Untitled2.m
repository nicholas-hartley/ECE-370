clear all;
close all;

t = [-5:0.001:5];
D11 = mod(11808942, 11);
D4 = mod(11808942, 4);

T = 2;
W = 1;
K = 23 + D11;

xt = FSWave(t, K, T, W);

plot(t, real(xt), '.-r', 'LineWidth', 2);

title("xt v.s. t and yt v.s. t");
xlabel('t');
ylabel('real(yt)');

hold on;

% % Experiment A
% yt = a_FSWave(t, K, T, W);
% 
% % Should be same as xt since X sub -k == X sub k
% plot(t, real(yt), '.-b', 'LineWidth', 1);

% % Experiment B
% yt = b_FSWave(t, K, T, W);
% 
% % Should be a shift left by t0
% plot(t, real(yt), '.-b', 'LineWidth', 1);

% % Experiment C
% yt = c_FSWave(t, K, T, W);
% 
% % Seems like applying odd and then scaling
% plot(t, real(yt), '.-b', 'LineWidth', 1);
 
% Experiment D
yt = d_FSWave(t, K, T, W);

%
plot(t, real(yt), '.-b', 'LineWidth', 1);

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

function [xt] = a_FSWave(t,K,T,W)
    syms x k
    expr = (1-4*x^2)*exp(-1i*2*pi*k*x/T);
    xt = zeros(1,length(t));
    for i = -K:K % Not assuming discrete time
        expr2 = subs(expr, k, -i);
        Xk = (1/T)*vpaintegral(expr2, x, -W/4, W/4);
        % Summing vectors of same size
        xt = xt + SUMCS(t, Xk, 4*pi()/T*i);
    end
end

function [xt] = b_FSWave(t,K,T,W)
    syms x k
    t0 = 0.65;
    expr = (1-4*x^2)*exp(-1i*2*pi*k*x/T);
    xt = zeros(1,length(t));
    for i = -K:K % Not assuming discrete time
        expr2 = subs(expr, k, i);
        Xk = (1/T)*vpaintegral(expr2, x, -W/4, W/4);
        Xk = Xk*exp(j*2*pi()*i*t0/T);
        % Summing vectors of same size
        xt = xt + SUMCS(t, Xk, 4*pi()/T*i);
    end
end

function [xt] = c_FSWave(t,K,T,W)
    syms x k
    expr = (1-4*x^2)*exp(-1i*2*pi*k*x/T);
    xt = zeros(1,length(t));
    for i = -K:K % Not assuming discrete time
        expr2 = subs(expr, k, i);
        Xk = (1/T)*vpaintegral(expr2, x, -W/4, W/4);
        Xk = Xk*j*i*2*pi()/T;
        % Summing vectors of same size
        xt = xt + SUMCS(t, Xk, 4*pi()/T*i);
    end
end

function [xt] = d_FSWave(t,K,T,W)
    syms x k
    expr = (1-4*x^2)*exp(-1i*2*pi*k*x/T);
    xt = zeros(1,length(t));
    for i = -K:K % Not assuming discrete time
        expr2 = subs(expr, k, i);
        if i > 0
            expr2 = subs(expr, k, K+1-i);
        elseif i < 0
            expr2 = subs(expr, k, -(K+1+i));
        else
            expr2 = subs(expr, k, i);
        end
        Xk = (1/T)*vpaintegral(expr2, x, -W/4, W/4);
        % Summing vectors of same size
        xt = xt + SUMCS(t, Xk, 4*pi()/T*i);
    end
end
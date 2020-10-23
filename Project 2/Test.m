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
        xs = xs + A(i)*exp(j*omega(i)/2*t);
    end
end

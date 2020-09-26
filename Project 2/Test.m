clear all;
close all;

T = .001;

n = [0:T:20];

y = 1 - cos(n);

dydt = ((1-cos(n+1) - (1-cos(n)))/T);

d2ydt2 = ((((1-cos(n+2) - (1-cos(n+1)))/T)-((1-cos(n+1) - (1-cos(n)))/T))/T);

plot(n, d2ydt2 + y)
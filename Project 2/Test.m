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
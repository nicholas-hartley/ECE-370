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

% Problem 1-e
for i = 1:n
   %i
   tic
   a = [0:0.1:10]
   T(i) = toc;
   
   %ii
   tic
   a = [0:0.1:10];
   Z(i) = toc;
   
   %iii
   b = [];
   tic
   for j = 0:100
       b(j+1) = j*0.1;
   end
   Y(i) = toc;
   
   %iv
   b(101) = 0;
   tic
   for j = 0:100
       b(j+1) = j*0.1;
   end
   X(i) = toc;
end

i = mean(T)
ii = mean(Z)
iii = mean(Y)
iv = mean(X)

%% Problem 2

% Problem 2-a
clear all;
close all;
t = linspace(0,1,10);

x = cos(2*pi*t + pi/4);

figure(1);
plot(t,x);

% Problem 2-b
clear all;
close all;

t = linspace(0,1,10);
x = cos(2*pi*t + pi/4);

figure(1);
plot(t,x, 'b', 'LineWidth', 2);

hold on;


t = linspace(0,1,20);
x = cos(2*pi*t + pi/4);

plot(t,x, 'Color', '#77AC30', 'LineWidth', 2);
%Plotted a darker green so it was visible

t = linspace(0,1,100);
x = cos(2*pi*t + pi/4);

plot(t,x, 'r--', 'LineWidth', 2);
%Made this line dashed so it less heavily covered green line

% Problem 2-c

% ---heaviside.m---
% function f = heaviside(t)
%     f = (t >= 0);
% end
% -----------------

y = heaviside([-1:0.2:1])
% Unit-step/heaviside function

% Problem 2-d
clear all;
close all;

t = [-1:1/100:1];

Y = heaviside(t);
G1 = heaviside(-t);
G2 = heaviside(t + 0.5);
G3 = heaviside(t - 0.5);
G4 = heaviside(-2*t + 0.5);

%OG Heaviside
subplot(5,1,1);
plot(t, Y, 'LineWidth', 2);
grid on; grid minor;
title('heaviside(t)');

%Negative Time
subplot(5,1,2);
plot(t, G1, 'LineWidth', 2);
grid on; grid minor;
title('heaviside(-t)');

%Shift left
subplot(5,1,3);
plot(t, G2, 'LineWidth', 2);
grid on; grid minor;
title('heaviside(t+0.5)');

%Shift right
subplot(5,1,4);
plot(t, G3, 'LineWidth', 2);
grid on; grid minor;
title('heaviside(t - 0.5)');

%More transforms
subplot(5,1,5);
plot(t, G4, 'LineWidth', 2);
grid on; grid minor;
title('heaviside(-2*t + 0.5)');

% Problem 2-e
t = linspace(0, 10); % period is where pi*t/5 = 0, 2*pi so (0,10)

% Assuming y(t) is a function of the instantaneous power so |y(t)|^2 can 
% be summed to find the energy over a time t

y = cos(pi*t/5);

E = abs(y).^2;

TotalE = sum(E)
Power = sum(E)/10 % Average for one period is total energy/period

% Problem 2-f
clear all;
close all;

n = linspace(-5,15,50) % Didn't specify # of points

tiledlayout(3,2)

%i
%(1,1)
stem(nexttile, n, heaviside(n));

%ii
%(1,2)
stem(nexttile, n, n.*heaviside(n));

%iii
%(2,1)
stem(nexttile, n, (0.8.^n).*heaviside(n));

%iv
%(2,2)
stem(nexttile, n, ((-0.8).^n).*heaviside(n));

%v
%(3,1)
stem(nexttile, n, sin(pi*n/2).*heaviside(n));

%vi
%(3,2)
stem(nexttile, n, 0.9.^n.*(sin(pi*n/4) + cos(pi*n/4)));

%% Problem 3

% Each subsection of 3 is done in the function with a the corresponding
% name (i.e. a is the function for problem 3.a)

clear all;
close all;

t = linspace(0,5,500);

a = cos(2*pi*t) + sin(5*pi*t);
b = cos(2*pi*t) + sin(2*t);
c = cos(2*t);
d = cos(2*pi*t/3);
e = cos(pi*t.^2/2);

subplot(5,1,1);
plot(t, a, 'LineWidth', 2);
grid on; grid minor;
title('cos(2*pi*t) + sin(5*pi*t)');

subplot(5,1,2);
plot(t, b, 'LineWidth', 2);
grid on; grid minor;
title('cos(2*pi*t) + sin(2*t)');

subplot(5,1,3);
plot(t, c, 'LineWidth', 2);
grid on; grid minor;
title('cos(2*n)');

subplot(5,1,4);
plot(t, d, 'LineWidth', 2);
grid on; grid minor;
title('cos(2*pi*n/3');

subplot(5,1,5);
plot(t, e, 'LineWidth', 2);
grid on; grid minor;
title('cos(pi*n^2/2)');

%% Problem 4

% Problem 4-a through d (only f0 changes in value)
clear all;
close all;

t = [0:1/8192:1];
f0 = 440;

x1 = cos(2*pi*f0*t);

subplot(5,1,1);
plot(t, x1, 'LineWidth', 2);
grid on; grid minor;
title('x1');

sound(x1);

% Problem 4-e and f
clear all;
close all;

t = [0:1/8192:1];
f0 = 440;
a = 2

x1 = cos(2*pi*f0*t);

x2 = exp(-(a^2+2)*t).*cos(2*pi*f0*t);

subplot(2,1,1);
plot(t, x1, 'LineWidth', 2);
grid on; grid minor;
title('x1');

soundsc(x1);

pause(3); %Time to hear difference between signals

subplot(2,1,2);
plot(t, x2, 'LineWidth', 2);
grid on; grid minor;
title('x2');

soundsc(x2);

% Problem 4-g and h
clear all;
close all;

t = [0:1/8192:1];
f0 = 510;
f1 = 12;
a = 2;

x1 = cos(2*pi*f0*t);

x2 = exp(-(a^2+2)*t).*cos(2*pi*f0*t);

x3 = cos(2*pi*f1*t).*cos(2*pi*f0*t);

x4 = 0.5*(cos(2*pi*t*(f1-f0)) + cos(2*pi*t*(f1+f0))); %Trig identity

subplot(3,1,1);
plot(t, x1, 'LineWidth', 2);
grid on; grid minor;
title('x1');

soundsc(x1);

pause(3); %Time to hear difference between signals

subplot(3,1,2);
plot(t, x2, 'LineWidth', 2);
grid on; grid minor;
title('x2');

soundsc(x2);

pause(3); %Time to hear difference between signals

subplot(3,1,3);
plot(t, x3, 'LineWidth', 2);
grid on; grid minor;
title('x3');

soundsc(x3);

pause(3);

soundsc(x4);
%% Problem 5

% Problem 5-a
clear all;
close all;

syms t;
f0 = 4e2;
phi = f0*t;
x1 = cos(2*pi*phi);
fins = matlabFunction(diff(phi)); % instantaneous frequency
fins() % call the function

% Problem 5-b
clear all;
close all;

syms t;
syms t0;
syms alpha;
phi = 0.5*alpha*t^2; % Times 1/2, so 2*pi*phi = pi*alpha*t^2
x4 = cos(2*pi*phi);
fins = matlabFunction(diff(phi)); % instantaneous frequency
case_t = fins(alpha, t) % All time case
case_0 = fins(alpha, 0) % instantaneous frequency evaluated at t = 0
case_t0 =fins(alpha, t0) % t = t0 case

% Problem 5-c
clear all;
close all;

t = [0:1/8192:1];
alpha = 1708;

x4 = cos(pi*alpha*t.^2);
soundsc(x4);

% Problem 5-d
clear all;
close all;

t = [0:1/8192:1];
alpha = 1708;

x4 = cos(pi*(alpha)*t.^2);
x4a1 = cos(pi*(alpha/2)*t.^2);
x4a2 = cos(pi*(alpha*2)*t.^2);
soundsc(x4);

pause(3);

soundsc(x4a1);

pause(3);

soundsc(x4a2);

% Problem 5-e
clear all;
close all;

syms t; 
time = [0:1/8192:1];

x5 = cos(2*pi*(exp(time)-100*time.^2+800*time));

soundsc(x5);

phi = exp(t)-100*t^2+800*t;

x5 = cos(2*pi*phi);

fins = matlabFunction(diff(phi)); % instantaneous frequency

t0 = fins(0)
t1 = fins(1)
t2 = fins(2)
%% Problem 6

clear all;
close all;
 
t = [0:1/8192:1];
alpha = 1708;

for phi = 0:pi/4:pi
    % phi % Uncomment to show the current phi for the loop
    x = cos(2*pi*alpha*t + phi);

    soundsc(x);

    pause(1);
end
%% Problem 7

%----notecreate.m----
% function [note] = notecreate(frq_no, stop)
% note = sin(2*pi* [1:stop]/8192 * (440*2.^((frq_no-1)/12)));
% end 
%--------------------

clear all;
close all;
 
Notename = {'A' 'A#' 'B' 'C' 'C#' 'D' 'D#' 'E' 'F' 'F#' 'G' 'G#'};
song = {'B' 'A' 'G' 'A' 'B' 'B' 'B' 'B' 'A' 'A' 'A' 'A' 'B' 'B' 'B' 'B' 'B' 'A' 'G' 'A' 'B' 'B' 'B' 'B' 'A' 'A' 'B' 'A' 'G' 'G'};

songidx(length(song)) = 0;

for k1 = 1:length(song) 
    idx = strcmp(song(k1), Notename); % Makes mask array of Notename for note
    songidx(k1) = find(idx); % Places index value into songidk to choose pitch
end
stop = 0.3*8192;
songnote = [ ];

for k1 = 1:length(songidx)
    songnote = [songnote; [notecreate(songidx(k1),stop) zeros(1,75)]'];
end
soundsc(songnote, 8192)

%plot(songnote)
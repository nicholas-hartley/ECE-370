%Nicholas Hartley

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

%soundsc(x1);

pause(3); %Time to hear difference between signals

subplot(3,1,2);
plot(t, x2, 'LineWidth', 2);
grid on; grid minor;
title('x2');

%soundsc(x2);

pause(3); %Time to hear difference between signals

subplot(3,1,3);
plot(t, x3, 'LineWidth', 2);
grid on; grid minor;
title('x3');

%soundsc(x3);

pause(3);

%soundsc(x4);
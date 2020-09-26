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

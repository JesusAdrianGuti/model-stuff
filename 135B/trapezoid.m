% this code use trapezoid method to approximate the integral of e^(-x^2) on
% the interval [0,1]
a = 0;
b = 1;
N = 58;
h = (b-a)/N;


f = @(x) exp(-x^2);

T = 0;

for i = 2:N+1
    T = T + (f(x(i-1))+f(x(i)))/2*h;
end

disp(T);
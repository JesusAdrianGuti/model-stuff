% simpson's 1/3 rule
a = 0;
b = 1;
f=@(x) exp(-x^2);

N = 100;
h = (b-a)/N;

x = a:h:b;

S = 0;

for i = 1:2:N-1
    S = S + h/3*( f(x(i)) + 4*f(x(i+1)) + f(x(i+2)) );
end

disp(S);
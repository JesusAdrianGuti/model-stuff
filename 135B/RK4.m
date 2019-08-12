% code RK4 method for u' = f(t,u)

f =@(t,u) 2*u;

t0 = 0;
tfinal = 10;
N = 10;
h = (tfinal-t0)/N;

t = t0:h:tfinal;
u = zeros(N+1,1);
u(1) = 0.5;

for i = 2:N+1
    b1 = 1/6;
    b2 = 2/6;
    b3 = 2/6;
    b4 = 1/6;
    k1 = f(t(i-1),u(i-1));
    k2 = f(t(i-1)+1/2*h,u(i-1)+h*1/2*k1);
    k3 = f(t(i-1)+1/2*h,u(i-1)+h*1/2*k2);
    k4 = f(t(i-1)+1*h, u(i-1)+h*1*k3);
    u(i) = u(i-1) + h*( b1*k1 + b2*k2 + b3*k3 + b4*k4);
end

uexact =@(t) u(1)*exp(2*t);
uexact_vec = uexact(t);

figure(1);
plot(t,u,'*-',t,uexact_vec,'o-')


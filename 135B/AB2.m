% code 2-step Adam-Bashforth method
f =@(t,u) 2*u;

t0 = 0;
tfinal = 1;
N = 100;
h = (tfinal-t0)/N;

t = t0:h:tfinal;
u = zeros(N+1,1);
u(1) = 0.5;

% first update by Euler to get enough previous steps
u(2) = u(1) + h*f(t(1),u(1));

for i = 3:N+1
    
    u(i) = u(i-1) + 3/2*h*f(t(i-1),u(i-1)) - 1/2*h*f(t(i-2),u(i-2));
    
end

uexact =@(t) u(1)*exp(2*t);
uexact_vec = uexact(t);

figure(1);
plot(t,u,'*-',t,uexact_vec,'o-')
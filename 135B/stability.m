% ODE: u' = -10u (stable)

a=0;
b=10;

u0 = 1;


f=@(t,u)-10*u;

h=1/2;
N=(b-a)/h;
T=a:h:b;
Ueuler=zeros(N+1,1);
Ueuler(1) = u0;
Ube=Ueuler;

% Euler method
for i = 2:N+1
    Ueuler(i) = Ueuler(i-1) + h*f(i*h,Ueuler(i-1));
end


% backward Euler method
for i = 2:N+1
    Ube(i) = Ube(i-1)/(1+10*h);
end

Uexact = u0*exp(-10*T);
figure(1);
subplot(1,2,1)
plot(T,Ueuler,'b*-',T,Uexact,'r');
title('Euler method')
subplot(1,2,2)
plot(T,Ube,'b*-',T,Uexact,'r');
title('Backward Euler method');

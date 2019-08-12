% Euler's method to solve IVP

t0 = 0;
tend = 1;
h = 0.001;

f = @(t,u) u;

u0 = 1;
u = u0;

t = 0;

for i = 1:1000
    
    upre = u;
    u = upre + h*f(t,upre);
    t = t + h;
   
%     disp([t,u])
    
end

ufun =@(t) exp(t);

disp([ufun(t) u])

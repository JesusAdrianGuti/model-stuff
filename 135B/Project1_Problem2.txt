% Math 135B Project 1
% Jesus Adrian Gutierrez
% Problem # 2
% Section 7-8pm

% Solving Initial Value Problem
% u(t) = ?4u(t) + 1
% u(0) = 0.5

% Predictor-Corrector Adams-Moulton VS. Runge-Kutta 4th order
format long

% ------------------ Using time step h = 0.5; ----------------------------

f=@(t,u) -4*u+1;

tini = 0.0;
tfin = 5.0;
h = 0.5;
n = (tfin - tini) / h;
k = 0;

t = tini : h : tfin;
u = zeros( n+1, 1);


% Find Exact Solution in order to compare our results
Exact = zeros(n+1 , 1);
for i = 1 : n+1
    Exact(i) = (1/4) + (1/4) * exp(-4 * t(i));
end


% Runge Kutta 4th Order Method
u(1) = 0.5;
for i = 1 : 3
    
    k1=h*f( t(i)         , u(i) );
    k2=h*f( t(i) + (h/2) , u(i)+(k1/2) );
    k3=h*f( t(i) + (h/2) , u(i)+(k2/2) );
    k4=h*f( t(i) + (h)   , u(i)+(k3) );

    u(i+1)= u(i) + (1/6) * (k1 + 2*k2 + 2*k3 +k4);
end


% Adams-Moulton Method 
L = zeros(n+1,1);
for i=5:n+1

L(i) = u(i-1)+(h/24)*(55*f(t(i-1),u(i-1))-59*f(t(i-2),u(i-2))+37*f(t(i-3),u(i-3))-9*f(t(i-4),u(i-4)));
u(i) = u(i-1)+(h/24)*(9*f(t(i),L(i))+19*f(t(i-1),u(i-1))-5*f(t(i-2),u(i-2))+f(t(i-3),u(i-3)));

end

figure(27)
% subplot(2,2,1)
plot(t,Exact,'b', t , u, '-*')
title('when h = 0.5')

% ------------------ Using time step h = 0.1; ----------------------------

f=@(t,u) -4*u+1;

tini = 0.0;
tfin = 5.0;
h = 0.1;  % <--------------- change to h = 0.1
n = (tfin - tini) / h;
k = 0;

t = tini : h : tfin;
u = zeros( n+1, 1);


% Find Exact Solution in order to compare our results
Exact = zeros(n+1 , 1);
for i = 1 : n+1
    Exact(i) = (1/4) + (1/4) * exp(-4 * t(i));
end


% Runge Kutta 4th Order Method
u(1) = 0.5;
for i = 1 : 3
    
    k1=h*f( t(i)         , u(i) );
    k2=h*f( t(i) + (h/2) , u(i)+(k1/2) );
    k3=h*f( t(i) + (h/2) , u(i)+(k2/2) );
    k4=h*f( t(i) + (h)   , u(i)+(k3) );

    u(i+1)= u(i) + (1/6) * (k1 + 2*k2 + 2*k3 +k4);
end


% Adams-Moulton Method 
L = zeros(n+1,1);
for i=5:n+1

L(i) = u(i-1)+(h/24)*(55*f(t(i-1),u(i-1))-59*f(t(i-2),u(i-2))+37*f(t(i-3),u(i-3))-9*f(t(i-4),u(i-4)));
u(i) = u(i-1)+(h/24)*(9*f(t(i),L(i))+19*f(t(i-1),u(i-1))-5*f(t(i-2),u(i-2))+f(t(i-3),u(i-3)));

end

figure(28)
%subplot(2,2,2)
plot(t,Exact,'k', t , u, '-*')
title('when h = 0.1')

% ------------------ Using time step h = 0.01; ----------------------------

f=@(t,u) -4*u+1;

tini = 0.0;
tfin = 5.0;
h = 0.01;
n = (tfin - tini) / h;
k = 0;

t = tini : h : tfin;
u = zeros( n+1, 1);


% Find Exact Solution in order to compare our results
Exact = zeros(n+1 , 1);
for i = 1 : n+1
    Exact(i) = (1/4) + (1/4) * exp(-4 * t(i));
end


% Runge Kutta 4th Order Method
u(1) = 0.5;
for i = 1 : 3
    
    k1=h*f( t(i)         , u(i) );
    k2=h*f( t(i) + (h/2) , u(i)+(k1/2) );
    k3=h*f( t(i) + (h/2) , u(i)+(k2/2) );
    k4=h*f( t(i) + (h)   , u(i)+(k3) );

    u(i+1)= u(i) + (1/6) * (k1 + 2*k2 + 2*k3 +k4);
end


% Adams-Moulton Method 
L = zeros(n+1,1);
for i=5:n+1

L(i) = u(i-1)+(h/24)*(55*f(t(i-1),u(i-1))-59*f(t(i-2),u(i-2))+37*f(t(i-3),u(i-3))-9*f(t(i-4),u(i-4)));
u(i) = u(i-1)+(h/24)*(9*f(t(i),L(i))+19*f(t(i-1),u(i-1))-5*f(t(i-2),u(i-2))+f(t(i-3),u(i-3)));

end

figure(29)
% subplot(2,2,3)
plot(t,Exact,'r', t , u, '-*')
title('when h = 0.01')


% ------------------ Using time step h = 0.001; ----------------------------


f=@(t,u) -4*u+1;

tini = 0.0;
tfin = 5.0;
h = 0.001;
n = (tfin - tini) / h;
k = 0;

t = tini : h : tfin;
u = zeros( n+1, 1);


% Find Exact Solution in order to compare our results
Exact = zeros(n+1 , 1);
for i = 1 : n+1
    Exact(i) = (1/4) + (1/4) * exp(-4 * t(i));
end


% Runge Kutta 4th Order Method
u(1) = 0.5;
for i = 1 : 3
    
    k1=h*f( t(i)         , u(i) );
    k2=h*f( t(i) + (h/2) , u(i)+(k1/2) );
    k3=h*f( t(i) + (h/2) , u(i)+(k2/2) );
    k4=h*f( t(i) + (h)   , u(i)+(k3) );

    u(i+1)= u(i) + (1/6) * (k1 + 2*k2 + 2*k3 +k4);
end


% Adams-Moulton Method 
L = zeros(n+1,1);
for i=5:n+1

L(i) = u(i-1)+(h/24)*(55*f(t(i-1),u(i-1))-59*f(t(i-2),u(i-2))+37*f(t(i-3),u(i-3))-9*f(t(i-4),u(i-4)));
u(i) = u(i-1)+(h/24)*(9*f(t(i),L(i))+19*f(t(i-1),u(i-1))-5*f(t(i-2),u(i-2))+f(t(i-3),u(i-3)));

end

figure(26)
% subplot(2,2,4)
plot(t,Exact,'r', t , u, '-*')
title('when h = 0.001')

% --------------------------------------------------------------------------




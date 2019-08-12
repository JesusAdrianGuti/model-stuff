% Math 135B Project 1
% Jesus Adrian Gutierrez
% Problem # 1
% Section 7-8pm

% To see more digits from each solution
format long

disp(' Numerically integrate f(x) from (0,1) ')
disp(' f(x) = x * e^(-x) ')
disp(' ')

f = @(t) t * exp(-t);
% Initial time t
tini = 0;
% Final time t
tfin = 1;
% Step size
h = 0.05;

N = ( tfin - tini )/ h;
t = tini: h : tfin;


% Euler Method for Numerical Integration
Euler = 0.0;
 for i = 1:N
     Euler = Euler + h * f(t(i));
 end
 disp('Euler = ');
 disp(Euler);
 
 
% Trapezoid Method
 Trapezoid = 0.0;
 for i = 2 : N+1
    Trapezoid = Trapezoid + ( f(t(i-1)) + f(t(i)) )/ 2 * h;
 end
disp('Trapezoid = ');
disp(Trapezoid);


% Composite Simpson's 1/3 Rule
Simpson =0.0;
simp1 = 0.0;
simp2 = 0.0;
    %find previous points
for i = 1 : N/2
    t(i) = tini + h * (2 * (i-1) );
    simp1 = simp1 + f(t(i));
end

for i = 1: (N-2)/2
    t(i) = tfin + 2 * h * i;
    simp2 = simp2 + f(t(i));
end
Simpson = (h/3) *( f(tini) + f(tfin) +  4*simp1 + 2*simp2);  
disp('Simpson = ');
disp(Simpson);


%  Gaussian Composite Quadrature Rule
Quadrature = 0.0;
t = tini: h : tfin;
node = [ -sqrt(1/3) sqrt(1/3)];
w = [1 1];
for i = 1 : N
        for k = 1:2
            Quadrature = Quadrature + (h/2) * w(k) * f( h/2 * (node(k)+1) + t(i));
        end
end
disp('Quadrature = ');
disp(Quadrature);

% Finally we find the Exact Solution in order to compare our results

Exact = 0.0;
Exact = integral(f,0,1,'ArrayValued',1);
disp('Exact Solution is ')
disp(Exact);





 
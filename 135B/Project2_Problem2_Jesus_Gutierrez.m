% Math 135B Project @
% Jesus Adrian Gutierrez
% Problem # 2
% Section 7-8pm


% Least square method. Finding coeffieicnts of for chebyshevs polynomial
% approximation
%  c0*T0(z) + c1T1(z) + c2T2(z)

% Data Points
x = [-2 -1 0 1 2];
y = [ 2  1 1 1 2];

% interval
a = min(x);
b = max(x);

% redefining our variable
z = (2*x-a-b)/(b-a);
n = 2;
[~,m]=size(z); 


for k = 1:m
    % define/construct chebyshev Polynomials
T(1,k) = 1;
T(2,k) = z(k);

for j = 3:n+1
T(j,k) = 2*z(k)*T(j-1,k)-T(j-2,k);
end

end

% X' * X * C = X' * Y
% 
A = T * T';
b = T * y';
c = A/b';

disp(" Chebyshev Polynomial Approximation ");
disp(" c0*T0(z) + c1T1(z) + c2T2(z)  "); 
disp(" coeffiencients are ");
disp(c);

% Math 135B Project @
% Jesus Adrian Gutierrez
% Problem # 1
% Section 7-8pm

% To see more digits from each solution
format long


% Using the Jacobi, Gauss-Seidel, SOR with w = 1.1 and W = 2 iterative methods, write
% and execute a computer program to solve the following linear system to four decimal places
% (rounded) of accuracy:

% Ax = b

A = [ 7  1 -1  2;
      1  8  0 -2;
     -1  0  4 -1;
      2 -2 -1  6; ];
       
b = [3 -5 4 -3]';


%================= Jacobi Method ======================

kmax = 100;
delta = 10^-10;
epsilon = 1e-4;

n = size(A,1);
x = zeros(n,1);
y = x; % previous x vector


for k = 1:kmax
    y = x;
    for i = 1:n
       
       % diag = A(i,i);
       % if ( abs(diag) < delta )
       %     disp("Diagonal element too small");
       %end
       
       x(i) = -(sum(A(i,:)'.*y)-A(i,i)*y(i))/A(i,i)+b(i)/A(i,i); 
    
    end
    
    if ( norm(x - y) < epsilon )
        disp(" ")
        disp(" Jacobi Method");
        disp(" Number of Iteraions "); 
        disp(k);
        disp(" Solution ");
        disp(x);
        break
    end
end


%================= Gauss-Seidel Method =========================

x = zeros(n,1);
for k = 1:kmax
    
    for i = 1:n
        
      x(i) = -(sum(A(i,:)'.*x)-A(i,i)*x(i))/A(i,i)+b(i)/A(i,i); 
    
    end
    
    if ( norm(x - y) < epsilon )
        disp(" ")
        disp(" Gauss-Seidel method ");
        disp(" Number of Iteraions "); 
        disp(k);
        disp(" Solution ");
        disp(x);
        break
       
    end
    
end


%=================== SOR with w = 1.1 and W = 2 ==================

w = 1.1;
x= zeros(n,1);
k = 1;

for j = 1:kmax
    i = 1;
    for i = 1:n
        
      x(i) = w*(-(sum(A(i,:)'.*x)-A(i,i)*x(i))/A(i,i)+b(i)/A(i,i)) + (1-w)*x(i);
  
    end
    
     if ( norm(x - y) < epsilon *10 )
        disp(" ")
        disp(" SOR Method, w =1.1 ");
        disp(" Number of Iteraions "); 
        disp(j);
        disp(" Solution ");
        disp(x);
        break
    end
    
end


w = 2.0;
x= zeros(n,1);

for j = 1:kmax
    
    for i = 1:n
        
      x(i) = w*(-(sum(A(i,:)'.*x)-A(i,i)*x(i))/A(i,i)+b(i)/A(i,i)) + (1-w)*x(i);
  
    end
    
     if ( norm(x - y) < epsilon )
        w = 1.1;
x= zeros(n,1);
k = 1;

for j = 1:kmax
    i = 1;
    for i = 1:n
        
      x(i) = w*(-(sum(A(i,:)'.*x)-A(i,i)*x(i))/A(i,i)+b(i)/A(i,i)) + (1-w)*x(i);
  
    end
    
     if ( norm(x - y) < epsilon *10 )
        disp(" ")
        disp(" SOR Method, w =1.1 ");
        disp(" Number of Iteraions "); 
        disp(j);
        disp(" Solution ");
        disp(x);
        break
    end
    
end

        break
    end
    
end

if ( j == kmax)
    disp(" ")
    disp(" SOR Method, w = 2.0 ");
    disp(" Number of Iteraions "); 
    disp(j)
    disp(" Solution ");
    disp(x);
end










    
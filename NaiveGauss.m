function [x] = NaiveGauss(n,a,b)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

x=zeros(1,n);

for k = 1:n-1
    for i = k+1:n
        xmult = a(i,k) / a(k,k);
        a(i,k) = xmult;
        
        % update every entry of Row
        for j = k+1:n
            a(i,j) = ( a(i,j) - (xmult * a(k,j)) );
        end
        % update b vector same as matrix
        b(i) = b(i) - (xmult * b(k));
    end
end

% ===================  back substitutions  ========================

% the last equation is a simple in the form, coefficient * x = num
x(n) = b(n)/a(n,n);

% using the know Xi we will use that to get Xi-1 and so on
for i = n-1:-1:1
    sum = b(i);
    for j = i+1:n
        sum = sum - a(i,j)*x(j);
    end
    x(i) = sum/a(i,i)
end
    
    
    
    
    
    
    
end
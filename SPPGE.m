function [x,list] = math135_project1_Gutierrez_Jesus(size,matrixA,b)
% Ax = b 
% input size, n size of square matrix, matrix A, b 
% Scaled partial pivoting gaussian elimination
% the function solves systems of linear equations Ax = b
% using the SPPGE method. Solution is vector x and order list l

% initialize 
x = zeros(1,size);   % the solution vector
list = 1:1:size;        % list keeps track of the order of solution
s = zeros(1,size);   


% Define the S vector
for i = 1:size
    s(i) = max( abs( matrixA(i,:) ) );   % take the max value from row i of  matrix A
end


% Define ratio vector
for k = 1:size
    column = matrixA(:,k);
    for  i = 1:size
        ratio(i) = abs( column(i) ) / s(i) ;
    end
    [rmax rindex] = max(ratio);
    
    % update L to make rmax the first index
    list(k) = rindex;
    list(rindex) = k;
    
    for 
    
    
    
end
rmax = 
end

end


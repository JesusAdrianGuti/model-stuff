% composite Gaussian quadrature rule
% store the standard nodes and weights for [-1 1]
Nd = [-sqrt(3/5) 0 sqrt(3/5)];
Wt = [5/9 8/9 5/9];

a = 0;
b = 1;

N = 10;
h = (b-a)/N;
x = a:h:b;

f=@(x) exp(-x^2);

G = 0;

for i = 1:N
    
    for j = 1:3
        
        G = G + h/2*(Wt(j)*f(h/2*(Nd(j)+1)+x(i)));
        
    end
    
end
display(G)


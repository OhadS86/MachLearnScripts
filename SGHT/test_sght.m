m = 200;  
n = 200*25; 

% for reproducibility
sd = 13;

randn('state', sd);
A = randn(m,n);       
x0 = randn(n,1);
y = A * x0 + randn(m, 1);

ind = (0 : 100: n);
s2 = int32((length(ind)-1)/10);
s1 = s2 * 5;



tic
[x, fval] = sghtISTA(A, y, ind, s1, s2);
toc;

tic
[x_w, fval_w] = sghtISTAWolfe(A, y, ind, s1, s2);
toc


fprintf ('diff = %.5f\n', norm(x - x_w) / norm(x));
fprintf ('obj: fval = %.5f, fval_w = %.5f\n', fval(end), fval_w(end));

fprintf ('nnz: x = %d, x_w = %d, \n', nnz(x), nnz(x_w));

g = 0; gw = 0; 
for i = 1 : length(ind) - 1
    idx = ind(i) + 1 : ind(i+1);
    if norm(x(idx)) > 0
        g = g + 1;
    end
    if norm(x_w(idx)) > 0
        gw = gw + 1;
    end
end

fprintf ('nng: x = %d, x_w = %d\n', g, gw);

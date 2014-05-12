function [ x, fval ] = sghtISTA( A, y, ind, s1, s2 )
% sghtISTA Summary of this function goes here
%   Detailed explanation goes here

[m, n] = size(A);
grp = length(ind) - 1;
x = zeros(n, 1);

maxIter = 1000;
fval = zeros(maxIter, 1);
tol = 1e-5;

%% Sanity check
assert(ind(grp+1) == n);
assert(ind(1) == 0);
assert(s1 <= n);
assert(s2 <= grp);

%%
fval(1) = 0.5 * norm(A * x - y, 2)^2;
x_last = x;
g_last = A' * (A * x - y);
L = 1;

while true
    x_cur = sght(x_last - g_last/L, ind, s1, s2);
    ax = A * x_cur;
    fval_cur = 0.5 * norm(ax - y, 2)^2;
    if (fval_cur < fval(1) - L/2 * norm(x_cur - x_last, 2)^2)
        break
    end
    L = L * 2;
%         fprintf ('\t\t L = %f\n', L);
end

fval(2) = fval_cur;


for itr = 3 : maxIter
    
    g_cur = A' * (ax - y);
    
    delta_x  = x_cur - x_last;
    delta_g = g_cur - g_last;
    
    if (norm(delta_x) == 0 || delta_x' * delta_g == 0)
        L = 1;
    else
        L = (delta_x' * delta_g) / (delta_x' * delta_x);
    end
    
%     if (rem(itr, 10) == 0)
%         fprintf ('\tIteration %d\n', itr);
%         fprintf ('\t\tL is initialized to %f\n', L);
%         fprintf ('\t\tobj is %f\n', fval(itr-1));
%     end
        
    while true
        x = sght(x_cur - g_cur/L, ind, s1, s2);
        ax = A * x;
        f = 0.5 * norm(ax - y, 2)^2;
        if (f < fval(itr-1) - L/2 * norm(x - x_cur, 2)^2)
            break
        end
        L = L * 2;
    end
    fval(itr) = f;
    x_last = x_cur;
    g_last = g_cur;
    x_cur = x;
    
    if (abs(f - fval(itr-1)) < tol * fval(itr-1) || norm(g_cur) < tol || f < tol)
        break
    end
end

fval = fval(1 : itr);
end


function [ x_admm ] = sglp_admm( v, s1, s2, ind, f_cmp )
%SGLP_ADMM Summary of this function goes here
%   Using ADMM to solve:
%    \min_{x}      1/2 \| x - v \|_2^2
%    subject to    \|x\|_1 <= s1
%                  \sum \|x_{G_i}\|_2 <= radius,
%
%  where the group information is contained in the "ind" array
%  Require SLEP package, the "eplb" function.


n = length(v);
grpNum = length(ind) - 1;
tol = 1e-8;

%% Temporary Arrays Only used for glLeastC
tmpU = zeros(grpNum, 1);
tmpW = zeros(grpNum, 1);

%% Main Procedure
if norm(v, 1) <= s1 && gnorm(v, ind) <= s2
    x_admm = v;
    return;
end
	
%% ADMM Initialization:

% Penalty
rho = 100;
rho_last = rho;

% multiplier: 
% lbd: n by 1 vector
% eta : n by 1 vector
lbd = ones(n, 1);
eta = ones(n, 1);
u = ones(n, 1);
w = ones(n, 1);
x = zeros(n, 1);
% exp = ones(n, 1);
iterNum = 0;

%% ADMM Iteration
while true

    %% fix lbd and nv, solve x
    % system of linear equation
    x_last = x;
    u_last = u;
    w_last = w;
%     lbd_last = lbd;
%     eta_last = eta;

    %% fix u and w, solve x
    x = 1/(1 + 2 * rho) * (v + rho * (u_last + w_last + lbd + eta));
    
    %% fix x and u, solve w: \ell_{2,1} ball projection
	[w, ~] = glLeastC(x - eta, s2, ind, tmpW, tmpU );

	%% fix x and w, solve u: \ell_1 ball  projection
    [u, ~, ~]=eplb(x - lbd, n, s1, 0);
    
    %% Update the multiplier
    lbd = lbd + u - x;
    eta = eta + w - x;
    
	%% Stopping Criterion
    iterNum = iterNum + 1;
    if mod(iterNum, 10000) == 0
		fprintf('Iteration %d:\n', iterNum);
		fprintf('Deviation x - u = %.6f, x - w= %.6f gnorm = %.3f, 1-norm = %.3f\n', norm(x-u), norm(x-w), gnorm(u, ind), norm(u,1));
		%fprintf('Last: u = %.3f, v=%.3f, x=%.3f, lbd=%.3f, nu = %.3f,  rho = %d\n\n\n', norm(u-u_last), norm(w-w_last), norm(x-x_last), norm(lbd-lbd_last), norm(nu-nu_last), rho);
    end
	
    f_admm = 0.5 * norm(x - v)^2;
    if norm(x, 1) <= s1 && gnorm(x, ind) <= s2 && f_admm < f_cmp + 0.1
        break;
    end
    if norm(x - u) < tol && norm(x - w) < tol
        break
    elseif norm(x-x_last) < tol * norm(x_last) && norm(u-u_last) ...
			< tol * norm(u_last) && norm(w-w_last)...
			< tol * norm(w_last)
        rho_last = rho;
        rho = rho * 2;
        lbd = lbd * 2;
        eta = eta * 2;
    end
end
x_admm = x;
end


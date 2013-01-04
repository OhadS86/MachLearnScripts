function [ x, fval ] = sglLeastC_admm( A, y, s1, s2, opts )
%SPG_admm Summary of this function goes here
%   Using the ADMM to solve the Sparse Group Projection (SPG) problem:
%
%   minimize     \| Ax - y \|
%   subject to   \| x \|_1 <= s1
%                sum \| x \|_{G_i} <= s2


ind = opts.ind;
grpNum = length(ind) - 1;
[m n] = size(A);
%% ADMM Initialization:

%% Temporary Arrays Only used for glLeastC
tmpU = zeros(grpNum, 1);
tmpW = zeros(grpNum, 1);

% Penalty
rho = 10000;
[L U] = factor(A, 2 * rho);
rho_last = rho;

% multiplier: 
% lbd: n by 1 vector
% nv : n by 1 vector
lbd = ones(n, 1);
nu = ones(n, 1);
u = ones(n, 1);
v = ones(n, 1);
x = zeros(n, 1);
% exp = ones(n, 1);
iterNum = 0;

%% ADMM Iteration
while true

    %% fix lbd and nv, solve x
    % system of linear equation
    x_last = x;
    u_last = u;
    v_last = v;
    lbd_last = lbd;
    nu_last = nu;
%     exp_last = exp;

    %% fix u and v, solve x
% tm = cputime;
	if rho ~= rho_last
		[L U] = factor(A, 2 * rho);
	end
		exp = A'*y + rho*(u + lbd + v + nu);
    
	if( m >= n )    % if skinny
		x = U \ (L \ q);
	else            % if fat
		x = exp/(2*rho) - (A'*(U \ ( L \ (A*exp) )))/(2*rho)^2;
	end
%     x = (A'*A + 2 * rho * eye(n)) \ (exp);
% fprintf('Step 1: %.5f\n', cputime - tm);
    
	[v, ~] = glLeastC(x-nu, s2, ind, tmpU, tmpW);
	tt = gnorm (v, ind);
% fprintf('Step 2: %.5f\n', cputime - tm);
   

    %% fix x and v, solve u
    % \ell_1 ball  projection, using function 'eplb' from SLEP.
% tm = cputime;
	[u, lambda, zf_step]=eplb(x - lbd, n, s1, 0);
% fprintf('Step 3: %.5f\n', cputime - tm);
    
	%% Update the multiplier
	lbd = lbd + u - x;
	nu = nu + v - x;
    
    
    
	%% Stopping Criterion
	iterNum = iterNum + 1;
	if(mod(iterNum, 10000) == 0)
		fprintf('Iteration %d:\n', iterNum);
		fprintf('Deviation x - u = %.6f, x-v= %.6f norm(2,1) = %.3f, norm(1) = %.3f\n', norm(x-u), norm(x-v), tt, norm(u,1));
		fprintf('Last: u = %.3f, v=%.3f, x=%.3f, lbd=%.3f, nu = %.3f,  rho = %d\n\n\n', norm(u-u_last), norm(v-v_last), norm(x-x_last), norm(lbd-lbd_last), norm(nu-nu_last), rho);
	end
	if norm(x - u) < 1e-5 & norm(x-v) < 1e-5
		break
	elseif norm(x-x_last) < opts.relTol & norm(u-u_last) < opts.relTol & ...
          	norm(v-v_last) < opts.relTol
		rho_last = rho;
		rho = rho * 2;
		lbd = lbd * 2;
		nu = nu * 2;
	end
end

x(abs(x)<1e-5) = 0;
fval = norm(A*x-y);
end


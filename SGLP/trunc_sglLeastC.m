function [ x, f_final_val ] = trunc_sglLeastC( A, y, s1, s2, tau, opts )
%TRUNC_SGLLEASTC Summary of this function goes here
%     \min_{x}      1/2 * \| Ax - y \|_2 
%	  subject to    \sum_j J_tau (|x_j|) <= s1
%					\sum_j J_tau (\|x_G_j\|_2) <= s2
    
%   J_tau (z) = min( |z|/tau, 1)

if (nargin <5)
    error('\n Inputs: A, y and s should be specified!\n');
elseif (nargin==5)
    opts=[];
end

% Get the size of the matrix A
[m,n]=size(A);

% Verify the length of y
if (length(y) ~=m)
    error('\n Check the length of y!\n');
end

% Verify the value of z
if (s1 <= 0 || s2 <= 0)
    error('\n s1 and s2 should be positive!\n');
end

if ~isfield(opts,'ind')
    error('\n In trunc_sglLeastC, .ind should be specified');
end

if ~isfield(opts,'x0DC')
    error('\n In trunc_sglLeastC, .x0DC should be specified');
end

if ~isfield(opts,'maxIterDC')
    error('\n In trunc_sglLeastC, .maxIterDC should be specified');
end

ind = opts.ind;
grpNum = length(ind) - 1;
maxIterDC = opts.maxIterDC;
x0 = opts.x0DC;
tol = opts.tolDC;

%% Initialization
fval = zeros(maxIterDC, 1);
opts.supp1 = (find(abs(x0) <= tau));
opts.suppG = [];
for i = 1 : grpNum
	if norm(x0(ind(i)+1: ind(i+1))) <= tau
		opts.suppG = [opts.suppG; i];
	end
end

opts.init = 1;

%% Iteration
iterNum = 0;
x_last = x0;
x_k = x0;
f_last = inf;
while iterNum < maxIterDC
	last_idx = iterNum;
	iterNum = iterNum + 1;
	if s1 < n - length(opts.supp1) || s2 < grpNum - length(opts.suppG)
		break;
    end
    opts.x0 = x_last;
	[ x_k, ~ ] = trunc_DC_sglLeastC( A, y, ...
		tau * (s1 - n + length(opts.supp1)), ...
		tau * (s2 - grpNum + length(opts.suppG)), opts );
	
	fval(iterNum) = 0.5 * norm(A*x_k - y)^2 + 0.5 * opts.rsL2 * norm(x_k)^2;
	if f_last - fval(iterNum)  < tol
		break;
	end
% 	fprintf('DC iteration %d, fval=%.6f\n', iterNum, fval(iterNum));
% 	if norm(x_k - x_last) / norm(x_k) < tol
% 		break;
% 	end
	opts.supp1 = [];
	opts.supp1 = (find(abs(x_k) <= tau));
	opts.suppG = [];
	for i = 1 : grpNum
		if norm(x_k(ind(i)+1: ind(i+1))) <= tau
			opts.suppG = [opts.suppG; i];
		end
	end
	x_last = x_k;
	f_last = fval(iterNum);
end
x = x_last;
x(abs(x) <= tau) = 0;
f_final_val = fval(1 : last_idx);
% Error('end here');
end


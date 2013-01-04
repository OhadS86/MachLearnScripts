function [ x_dykstra] = sglp_dykstra( v, s1, s2, ind, f_cmp )
%SGLP_DYKSTRA Summary of this function goes here
%   Using Dykstra Algorithm to solve:
%    \min_{x}      1/2 \| x - v \|_2^2
%    subject to    \|x\|_1 <= s1
%                  \sum \|x_{G_i}\|_2 <= radius,
%
%  where the group information is contained in the "ind" array
%  Require SLEP package, the "eplb" function.
%  
%  Reference:
%    "Proximal Splitting Methods in Signal Processing" 
%           -- P. L. Combettes and J-C Pesquet
%

n = length(v);
grpNum = length(ind) - 1;
tol = 1e-5;
%% Temporary Arrays Only used for glLeastC
tmpU = zeros(grpNum, 1);
tmpW = zeros(grpNum, 1);

%% Main Procedure
if norm(v, 1) <= s1 && gnorm(v, ind) <= s2
    x = v;
    return;
else
	x_last = v;
	p_last = zeros(n, 1);
	q_last = zeros(n, 1);

	maxIter = 5000;

	for iter = 1 : maxIter
		[y_last, ~] = glLeastC( x_last + p_last, s2, ind, tmpW, tmpU );
		p_new = x_last + p_last - y_last;
		[x_new, ~, ~]=eplb(y_last + q_last, n, s1, 0);
		q_new = y_last + q_last - x_new;
		
        f_dykstra = 0.5 * norm(x_new - v)^2;
        if norm(x_new, 1) <= s1 && gnorm(x_new, ind) <= s2 ...
                && f_dykstra < f_cmp + 0.1
            break;
        end
		if norm(x_new-x_last) < tol * norm(x_new) && norm(p_new-p_last) ...
				< tol * norm(p_last) && norm(q_new-q_last) ...
				< tol * norm(q_last)
			break;
		end
		
	end
	x_dykstra = x_new;
end


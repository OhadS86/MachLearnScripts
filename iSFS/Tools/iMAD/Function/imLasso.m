function [x, funVal] = imLasso (AA, y, W, lambda, opts)
%Core optimization procedure called by iSFSLS_helper. 
%Essentially a lasso problem with weighting scheme

global ind S pf hasPf

% Get the dimension of beta
pfN = sum(hasPf);

% Verify the value of lambda
if (lambda < 0)
	error('\n lambda should be nonnegative!\n');
end

% run sll_opts to set default values (flags)
opts=sll_opts(opts);


%% May add rsL2 later Modification Needed (4)
rsL2 = 0;

%% The main program

%% The Armijo Goldstein line search scheme + accelearted gradient descent


bFlag=0; % this flag tests whether the gradient step only changes a little
L = 1 + rsL2;
n = size (AA, 2);

% assign xp with x
x = zeros (n, 1);
xp=x;
xxp = zeros (n,1);

% ap and aa are used for computing the weight in forming search point
ap=0; aa=1;

tot = 2^S - 1;


%% Precompute the A_mask matrix
A = AA;
A (isnan(A)) = 0;
T = cell (tot, 1);
for mask = 1 : tot
	if (~hasPf (mask))
		continue
	end
	% Extract all data that CONTAINS mask in its profile
	smp_idx =  (bitand (pf, mask) == mask);
	coeff = zeros (1, n);	
	for i = 1 : S
		coeff(ind(i)+1:ind(i+1))= repmat(W(mask, i), 1, ind(i+1) - ind(i));
	end
	T{mask} = A(smp_idx, :) * diag (coeff);
end

for iterStep = 1 : opts.maxIter
	% --------------------------- step 1 ---------------------------
	% compute search point s based on xp and x (with bt)
	bt = (ap-1)/aa;    s = x + bt* xxp;
	
	% --------------------------- step 2 ---------------------------
	% line search for L and compute the new approximate solution x
	
	%% compute the gradient (g) at s, Modification needed (1)
	g = zeros (n, 1);
	
	for mask = 1 : tot
		% no such profile
		if (~hasPf (mask))
			continue
		end
		smp_idx =  (bitand (pf, mask) == mask);
		g = g + T{mask}' * (T{mask} * s - y(smp_idx)) / (pfN * sum(smp_idx));
	end
	
	% copy x to xp
	xp = x;
	
	while (1)
		% let s walk in a step in the antigradient of s to get v
		% and then do the l1-norm regularized projection
		v = s - g / L;
		
		% L1-norm regularized projection
		x= sign(v) .* max(abs(v)-lambda / L,0);
		
		v= x - s;  % the difference between the new approximate solution x and the search point s
		
		% computel_sum, r_sum for line search criterion
		r_sum=v'*v;
		l_sum = 0;
		for mask = 1 : tot
			% no such profile
			if (~hasPf (mask))
				continue
			end
			% Extract all data that CONTAINS mask in its profile
			smp_idx =  (bitand (pf, mask) == mask);
			l_sum = l_sum + norm (T{mask} * v) ^2 / (pfN * sum(smp_idx));
		end
		
		
		if (r_sum <=1e-20)
			bFlag=1; % this shows that, the gradient step makes little improvement
			break;
		end
		
		% line search condition: l_sum <= (L - rsL2) * ||v||_2^2,
		% Modification Needed (2)
		if(l_sum <= r_sum * (L-rsL2))
			break;
		else
			L=max(2*L, l_sum/r_sum + rsL2);
		end
	end
	
	ValueL(iterStep)=L;
	
	% --------------------------- step 3 ---------------------------
	% update aa and ap, and check whether converge
	ap = aa; aa = (1+ sqrt(4 * aa * aa + 1))/2;
	
	xxp = x - xp;
	
	%% Compute function value
	ret = 0; 
	for mask = 1 : tot
		% no such profile
		if (~hasPf (mask))
			continue
		end
		
		% Extract all data that CONTAINS mask in its profile
		smp_idx =  (bitand (pf, mask) == mask);
		ret = ret + norm (T{mask} * x - y(smp_idx))^2 / (pfN * sum(smp_idx));
	end
	funVal(iterStep) = 0.5 * ret + sum(abs(x)) * lambda;
	
	if (bFlag)
		% fprintf('\n The program terminates as the gradient step changes the solution very small.');
		break;
	end
	
	switch(opts.tFlag)
		case 0
			if iterStep>=2
				if (abs( funVal(iterStep) - funVal(iterStep-1) ) <= opts.tol)
					break;
				end
			end
		case 1
			if iterStep>=2
				if (abs( funVal(iterStep) - funVal(iterStep-1) ) <=...
						opts.tol* funVal(iterStep-1))
					break;
				end
			end
		case 2
			if ( funVal(iterStep)<= opts.tol)
				break;
			end
		case 3
			norm_xxp=sqrt(xxp'*xxp);
			if ( norm_xxp <=opts.tol)
				break;
			end
		case 4
			norm_xp=sqrt(xp'*xp);    norm_xxp=sqrt(xxp'*xxp);
			if ( norm_xxp <=opts.tol * max(norm_xp,1))
				break;
			end
		case 5
			if iterStep>=opts.maxIter
				break;
			end
	end
end

end



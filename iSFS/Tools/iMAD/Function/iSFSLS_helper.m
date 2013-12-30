function model = iSFSLS_helper (X_local, Y_local, ind_local, S_local, pf_local, hasPf_local, lambda, config)
%Helper function: will be called by iSFSLS
%	Take care of detailed parameter settings.

% Input:
% S:		Scalar.						Number of Sources
% tot:		Scalar = 2^S - 1.			Total number of all combinations of groups
% n:		Scalar.						Number of samples
% p:		Scalar.						Number of features
% X:		n by p matrix.				The original data matrix. Possibly come with missing entry / group.
% pf:		n by 1 vector.				Profile of each sample. pf[i] & (1<<(j-1) != 0 means group j exists in the ith sample
% hasPf:	tot by 1 boolean vector.	Indicate whether certain profile exist in the data.
% lambda:	Scalar						regularization parameter

% helper functions
% updateBeta ().						Update beta. Essentially a Lasso problem
% updateAlpha ().						Update Alpha. This should be down for each profile independently Possibly using LeastC/LogisticC in SLEP package.

% output:
%		model.Alpha:    tot by S matrix.            Identical to W except Alpha is indexed by mask
%		model.beta:		p by 1 vector.				The linear model learned.
%		model.ind
%		model.S
%		model.hasPf
%		model.feaIdx:	p by 1 boolean				feature index



global X Y ind S pf hasPf beta

X = X_local; Y = Y_local; ind = ind_local; S = S_local; pf = pf_local; hasPf = hasPf_local;
[~, p] = size (X);

if ~isfield (config, 'alpha') || ~isfield (config, 'beta')
	error ('\t In iSFSLS: the update method of alpha and beta must be specified\n');
end
% Add other possible initialization to beta later
beta = 0.001 * ones (p, 1);

if strcmp (config.alpha, 'uniform')
	tot = 2^S - 1;
	powOfTwo = 2 .^ (0 : S - 1); % powOfTwo(i) = 2 ^ (i - 1)
	Alpha = zeros (tot, S);
	
	for mask = 1 : tot
		if (~hasPf (mask))
			continue
		end
		num = sum (bitand (mask, powOfTwo) > 0);
		Alpha (mask, bitand (mask, powOfTwo) > 0) = 1/num;
	end
	
	if strcmp (config.beta, 'default')
		beta = updateBeta (X, Y, Alpha, lambda);
	elseif strcmp (config.method, 'BCD')
		updateBetaBCD (Alpha, lambda);
	end
else
 	maxIter = 30;
	beta_last = -1 * ones (p, 1);
	
	% Initialize Alpha
	tot = 2^S - 1;
	powOfTwo = 2 .^ (0 : S - 1); % powOfTwo(i) = 2 ^ (i - 1)
	Alpha = zeros (tot, S);
	for mask = 1 : tot
		if (~hasPf (mask))
			continue
		end
		num = sum (bitand (mask, powOfTwo) > 0);
		Alpha (mask, bitand (mask, powOfTwo) > 0) = 1/num;
	end
	
	% Calculate beta and alpha iteratively
	fv_last = inf;
	for iter = 1 : maxIter
		if strcmp (config.beta, 'default')
			[beta, fv] = updateBeta (X, Y, Alpha, lambda);
		elseif strcmp (config.method, 'BCD')
			updateBetaBCD (Alpha, lambda);
		end
		Alpha = updateAlpha (config.alpha);
		
		%fprintf ('\t Inner iteration %d\n', iter);
		if norm (beta - beta_last) < 1e-3 * norm(beta) ...
				|| abs (fv - fv_last) < 1e-5 * fv
		%	fprintf ('Converge after %d iterations\n', iter);
			break
		end
		beta_last = beta;
		fv_last = fv;
	end
end

model.beta = beta;
model.Alpha = Alpha;
model.ind = ind;
model.S = S;
model.hasPf = hasPf;
model.feaIdx = (beta ~= 0);

clear X Y ind S pf hasPf beta
end


%% Sub-procedure that updates W iteratively
function Alpha = updateAlpha (method)
%
%----------------------- LeastC options -----------------------
%change to default value later for a neat code
opts=[];

% Starting point
opts.init=2;            % starting from a zero point

% Termination criterion
opts.tFlag=5;          % run .maxIter iterations
opts.maxIter=100;      % maximum number of iterations

% Mormalization
opts.nFlag=0;         % without normalization

global X Y ind S pf hasPf beta

tot = 2^S - 1;
powOfTwo = 2 .^ (0 : S - 1); % powOfTwo(i) = 2 ^ (i - 1)
Alpha = zeros (tot, S);

for mask = 1 : tot
	if (~hasPf (mask))
		continue
	end
	
	% Extract all data that CONTAINS mask in its profile
	smp_idx =  (bitand (pf, mask) == mask);
	X_mask = X (smp_idx, :);
	Y_mask = Y (smp_idx);
	
	A = [];
	for i = 1 : S
		if bitand (mask, powOfTwo (i)) > 0
			idx = (ind(i) + 1 : ind(i + 1));
			A = [A  X_mask(:, idx) * beta(idx)];
		end
	end
	if strcmp (method, 'l1')
		[sol, ~, ~] = LeastC (A, Y_mask, 1, opts);
	elseif strcmp (method, 'l2')
		[sol, ~, ~] = LeastCL2 (A, Y_mask, 1, opts);
	else
		error ('Update of alpha must be ''l1'' or ''l2''\n');
	end
	% 	[sol, ~] = nnLeastC (A, Y_mask, 1, opts);
	Alpha (mask, bitand (mask, powOfTwo) > 0) = sol';
end
end

%% Sub-procedure that updates the entire beta
function [sol, f] = updateBeta (X, y, Alpha, lambda)
%----------------------- imLasso options -----------------------
%change to default value later for a neat code
opts=[];

% Starting point
opts.init=2;            % starting from a zero point

% Termination criterion
opts.tFlag=1;          % relative change in objective value
opts.maxIter=1000;      % maximum number of iterations
opts.tol = 1e-5;


% Mormalization
opts.nFlag=0;         % without normalization

[sol, fval] = imLasso (X, y, Alpha, lambda, opts);
f = fval(end);
end
%% Sub-procedure that updates beta iteratively
function [] = updateBetaBCD (Alpha, lambda)
%----------------------- imLassoBCD options -----------------------
%change to default value later for a neat code
opts=[];

% Starting point
opts.init=2;            % starting from a zero point

% Termination criterion
opts.tFlag=1;          % relative change in objective value
opts.maxIter=1000;      % maximum number of iterations
opts.tol = 1e-3;

% Mormalization
opts.nFlag=0;         % without normalization

global X S beta

% Update beta for each group iteratively
beta_last = -1 * ones (size(X, 2), 1);
for round = 1 : 20
	for src = 1 : S
		res = getResidual (src, Alpha);
		[~, ~] = imLassoBCD (src, X, res, Alpha, lambda, opts); % Update beta_src
	end
	if norm (beta - beta_last) < 1e-3 * norm(beta)
		fprintf ('BCD %d rounds', round);
		break
	end
	beta_last = beta;
end
end


%% Calculate the residual of each mask w.r.t src th source
function res = getResidual (src, Alpha)
global X Y ind S pf hasPf beta

powOfTwo = 2 .^ (0 : S - 1); % powOfTwo(i) = 2 ^ (i - 1)
tot = 2^S - 1;

res = cell(tot, 1);

for mask = 1 : tot
	% no such profile, or the profile does not contain src
	if (~hasPf (mask) || (mask & powOfTwo(src)) == 0)
		continue
	end
	
	% Extract all data that CONTAINS mask in its profile
	smp_idx =  (bitand (pf, mask) == mask);
	X_mask = X(smp_idx, :);
	Y_mask = Y(smp_idx);
	
	X_mask (isnan(X_mask)) = 0;
	
	ret = Y_mask;
	for i = 1 : S
		if i == src
			continue
		end
		fea_idx = (ind(i) + 1 : ind(i + 1));
		ret = ret - Alpha(mask, i) * X_mask (:, fea_idx) * beta (fea_idx);
	end
	res{mask} = ret;
end
end






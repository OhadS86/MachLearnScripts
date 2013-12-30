function [ model ] = iSFSLS ( A, y, ind, paras, config )
%The iSFS model with least squares loss function. 
%This procedure supports cross validation. By default it is 5-fold
%Reference:

m = sum(y==1);
tot = length (y);
paN = length (paras);

%% If multiple arguments are provided, then 5-fold cross validation is employed
% Otherwise use all the data to train
if paN > 1 
	foldN = 5;
	mm = min (sum(y==1), sum(y==-1));
	blkSize = (mm - (rem (mm, foldN))) / foldN;
	
	best_prec = -1;
	for paCnt = 1 : paN
		precSum = 0;
		for cnt = 1 : foldN
			if cnt ~= foldN
				TePosIdx = (cnt - 1) * blkSize + 1 : cnt * blkSize;
				TrPosIdx = setdiff ( (1:m), TePosIdx);
				
				TeNegIdx = m + ((cnt - 1) * blkSize + 1 : cnt * blkSize);
				TrNegIdx = setdiff ( (m+1:tot), TeNegIdx);
			else
				TePosIdx = (cnt - 1) * blkSize + 1 : m;
				TrPosIdx = setdiff ( (1:m), TePosIdx);
				
				TeNegIdx = m + (cnt - 1) * blkSize + 1 : tot;
				TrNegIdx = setdiff ( (m+1:tot), TeNegIdx);
			end
			cvTr = [A(TrPosIdx, :); A(TrNegIdx, :)]; cvTrY = [y(TrPosIdx); y(TrNegIdx)];
			cvTe = [A(TePosIdx, :); A(TeNegIdx, :)]; cvTeY = [y(TePosIdx); y(TeNegIdx)];
			
			stat = calStat (cvTr, ind);
			model = iSFSLS_helper (cvTr, cvTrY, ind, stat.S, stat.pf, stat.hasPf, paras(paCnt), config);
			[our_lab, ~] = predict ('RandomForest', cvTr, cvTrY, cvTe, cvTeY, model);
			precSum = precSum + sum (our_lab == cvTeY) / length (our_lab);
		end
		%fprintf ('testing parameter %.9f\n', paLst (paCnt));
		prec = precSum / foldN;
		if prec > best_prec
			best_prec = prec;
			best_pa = paras(paCnt);
		end
	end
else
	best_pa = paras(1);
end

stat = calStat (A, ind);
model = iSFSLS_helper (A, y, ind, stat.S, stat.pf, stat.hasPf, best_pa, config);
end


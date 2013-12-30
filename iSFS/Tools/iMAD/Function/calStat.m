function [ stat ] = calStat( A, ind )
%Compute the profile of given data. 
%Reference:

S = length(ind) - 1;
startPt = ind(1 : S) + 1;
pf = bi2de (~isnan(A(:, startPt)));
hasPf = false (2^S - 1, 1);
hasPf (pf) = true;

stat.pf = pf;
stat.hasPf = hasPf;
stat.S = S;
end


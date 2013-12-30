%% Validate new Model. No parallel setting in Matlab

clc, clear
repN = 10;
acc = zeros (repN, 1); auc = zeros (repN, 1);
acc_imsf = zeros (repN, 1); auc_imsf = zeros (repN, 1);
bslAcc = zeros (repN, 1);


paras = logspace(-5, -1, 5);
config = [];
config.alpha = 'l1';
config.beta = 'default';

doption = [];
doption.clustering = true;
doption.ratio = .5;

for rep = 1 : repN
	[trainX, testX, trainY, testY, ind, S, pf, hasPf] = processData ('ADNC', doption);
    
	model = iSFSLS (trainX, trainY, ind, 0.01, config);
	[our_lab, ~] = predict ('RandomForest', trainX, trainY, testX, testY, model);
	[ROCX, ROCY, ~, auc(rep)] = perfcurve(testY, our_lab, 1);
	acc(rep) = sum(our_lab==testY) / length (our_lab);
	
	fprintf ('Iteration %d\n', rep);
	fprintf ('iSFS: %d / %d = %.5f with %d features, auc=%.5f\n', sum(our_lab == testY), length (our_lab), acc(rep), sum(model.feaIdx), auc(rep));
end





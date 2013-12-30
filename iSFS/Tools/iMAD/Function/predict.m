function [p_lab, conf] = predict (method, trainX, trainY, testX, testY, model)
%General prediction fuction.

if strcmp (method, 'RandomForest')
	fidx = model.feaIdx;
	
	trainX (isnan(trainX)) = 0;
	testX (isnan(testX)) = 0;
	trainX = trainX(:, fidx');
	testX = testX(:, fidx');
	
	RFModel = classRF_train(trainX, trainY);
	p_lab= classRF_predict (testX, RFModel);
	conf = p_lab;
else
	error ('\t Error: prediction method is not supported.');
end


end

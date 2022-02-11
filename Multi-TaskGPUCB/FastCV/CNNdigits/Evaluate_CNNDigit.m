function [Y_eval]= Evaluate_CNNDigit(nextPt, nextTasksBatch, X,Y,indices,layers)

for i=1:max(indices)
test = (indices == i);
training = ~test;
XValidation = X(:,:,:,test);
XTraining = X(:,:,:,training);
YValidation = Y(test);
YTraining = Y(training);

options = trainingOptions('sgdm', ...
    'MaxEpochs',20, ...
    'ValidationData',{XValidation,YValidation}, ...
    'ValidationFrequency',3, ...
    'Verbose',false,...
    'ExecutionEnvironment','gpu',...
    'InitialLearnRate',10^nextPt(1),...
    'L2Regularization',10^nextPt(2),...
    'Momentum',nextPt(3));

rng(10)
[~, netinfo] = trainNetwork(XTraining,YTraining,layers,options);

accuracy(i,:) = max(netinfo.ValidationAccuracy);
end
Y_eval = accuracy(nextTasksBatch,:);
end
function [Y_eval]= Evaluate_Arcene(nextPt, nextTasksBatch, Data, Label,indices)

for i=1:max(indices)
test = (indices == i);
training = ~test;
SVMModel = fitcsvm(Data(training,:),Label(training,:),'Standardize',true,...
    'KernelFunction','gaussian','BoxConstraint', 10^nextPt(1),'KernelScale',10^nextPt(2));
[label_pred,~] = predict(SVMModel,Data(test,:));
accuracy(i,:) = sum(label_pred == Label(test,:))/numel(Label(test,:));
end
Y_eval = accuracy(nextTasksBatch,:);
end
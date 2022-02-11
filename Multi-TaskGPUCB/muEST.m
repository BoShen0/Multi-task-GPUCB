function Muval = muEST(x, GPmodel, numTasks, Task2Numeric)
% Copyright (c) 2020 Bo Shen
% this function is used to calculate the infer optimal solution for a given
% gaussian process model
  numDims = size(x, 2);  % the dimension of the input
  numData = size(x, 1);

if numTasks ==1
    [mu, ~,~] = predict(GPmodel,x);
    Muval = mu;
else
  mus = zeros(numData, numTasks);
  for i = 1:numTasks
 % get the mu and sigma for each task with the same x
    [mu, ~,~] = predict(GPmodel,[x repmat(Task2Numeric(i,:),numData,1)] );  
    mus(:,i) = mu;
  end
  Muval = sum(mus, 2);
end

end


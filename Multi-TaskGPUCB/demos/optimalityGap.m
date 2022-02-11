[min(result_covariance{7}),min(result_CGPUCB{7}),min(result_singleTask{7}),...
    min(result_simpleSum{7})]

for i = 1:n
initalValue(i) = mean(fun(params.optimalSolution,1:params.numTasks)- fun(XtrRaw(i,:),1:params.numTasks));
end
min(initalValue)
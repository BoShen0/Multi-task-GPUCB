%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% first round revision
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
Tanparams.noiselevel  = 0.2;
Tanparams.T = 50;
Tanparams.Repetition = 10;
Tanparams.Task2Numeric =[0.33 0.67; 0.67 0.33];
Tanparams.targetset = 0;
Tanparams.optimalSolution = [0 0];
Tanparams.bounds = repmat([-32.768 32.768],2,1);


fun = @(x,i) -Evaluate_ackley(x, i);

n =10; % should be even number
rng(0,'twister'); % For reproducibility
XtrRaw = boGetInitPts(Tanparams.bounds,n, []);
Tasktr = [repmat(1, n/2, 1); repmat(2, n/2, 1)];
for i = 1:n
Ytr(i,1) = fun(XtrRaw(i,:),Tasktr(i));
end

Xtr = XtrRaw;
tic
output = paperTan(Xtr, Tasktr, Ytr, fun, Tanparams);
toc
for iter =1:Tanparams.Repetition
record_Tan(iter,:) = output{iter}(:,end); 
end
plot(mean(record_simpleSum(:,1:50)))
hold on;
plot(mean(cummin(record_Tan,2)))

%% prepare data for simple regret
for iter = 1:10
record_simpleSum(iter,:) = summary_simpleSum{iter}{9}; 
record_covariance(iter,:) = summary_covariance{iter}{9}; 
record_CGPUCB(iter,:) = summary_CGPUCB{iter}{9}; 
record_singleTask(iter,:) = repelem(summary_singleTask{iter}{9},2); 
end

DATA =[mean(record_simpleSum);  mean(record_covariance); ...
    mean(record_CGPUCB)+1;mean(record_singleTask)]';

DATAStd =[std(record_simpleSum);  std(record_covariance); ...
    std(record_CGPUCB);std(record_singleTask)]';

DATA = [DATA(1:50,:) mean(cummin(record_Tan,2))'];
DATAStd = [DATAStd(1:50,:) std(cummin(record_Tan,2))'];

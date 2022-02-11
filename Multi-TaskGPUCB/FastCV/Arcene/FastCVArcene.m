% This experiment is for fast validation method
%% prepare data for fast cross validation
NumFold = 5;
rng(10)
indices = crossvalind('Kfold',200,NumFold);
Data = Data_feat(1:200,:);
Label = Data_label(1:200);


% tStart = tic; 
% SVMModel = fitcsvm(Data_feat(1:160,:),Data_label(1:160,:),'OptimizeHyperparameters','all' );
% tEnd = toc(tStart)
% 
% [label,score] = predict(SVMModel,Data_feat(161:200,:));
% accuracy = sum(label== Data_label(161:200,:))/numel(Data_label(161:200,:))
%% for two tasks BO with same optimal solution
% params.numTasks = NumFold;
% params.fmincon.gridSize = 10000;
% params.bounds = repmat([-3 3],2,1);
% params.BetaUpdate= 'simple';
% % params.BetaUpdate= 'complex';
% params.prob= 0.2; 
% params.NumSearch = 1;
% % params.Task2Numeric = [];
% params.KernelFun = 'ardmatern52';
% params.T = 100;
% params.BatchSize = 1;
% params.Task2Numeric = [0 1 1 1 1;1 0 1 1 1;1 1 0 1 1; 1 1 1 0 1; 1 1 1 1 0];
% % params.Task2Numeric = [1 0 0;0 1 0;0 0 1];
% params.Strategy = 'MT-GPUCB';
% % params.Strategy = 'Random';
% % params.optimalSolution = [0 0];
% params.NewAcq = 'SimpleSum';
% params.MultiEvaluation = 'On';
% n =20; % should be even number

params.numTasks = NumFold;
params.fmincon.gridSize = 10000;
params.bounds = repmat([-3 3],2,1);
params.BetaUpdate= 'simple';
% params.BetaUpdate= 'complex';
params.prob= 0.2; 
params.NumSearch = 1;
% params.Task2Numeric = [];
params.KernelFun = 'ardmatern52';
params.T = 100;
params.BatchSize = 1;
params.Task2Numeric = [0 1 1 1 1;1 0 1 1 1;1 1 0 1 1; 1 1 1 0 1; 1 1 1 1 0];
% params.Task2Numeric = [1 0 0;0 1 0;0 0 1];
params.Strategy = 'MT-GPUCB';
% params.Strategy = 'Random';
% params.optimalSolution = [0 0];
params.NewAcq = 'SimpleSum';
params.MultiEvaluation = 'On';

fun = @(x,i) Evaluate_Arcene(x, i, Data, Label,indices);
% fun([0 2],[ 1 2])
n =20; % should be even number
rng(0,'twister'); % For reproducibility
XtrRaw = boGetInitPts(params.bounds,n, []);
Tasktr = [repmat(1, n/5, 1); repmat(2, n/5, 1);repmat(3, n/5, 1);repmat(4, n/5, 1);repmat(5, n/5, 1)];
for i = 1:n
Ytr(i,1) = fun(XtrRaw(i,:),Tasktr(i));
end

% for simple sum 
result_simpleSum = myBO(XtrRaw,Tasktr,Ytr,fun, params);

% for covariance 
params.NewAcq = 'Covariance';
result_covariance = myBO(XtrRaw,Tasktr,Ytr,fun, params);

% GPUCB
params.Strategy = 'CGP-UCB';
params.NewAcq = 'CGP-UCB';
result_CGPUCB = myBO(XtrRaw,Tasktr,Ytr,fun, params);

% single task
params.Strategy = 'MT-GPUCB';
params.NewAcq = 'Covariance';
params.BatchSize = 5;
params.T = 20;
result_singleTask = myBO(XtrRaw,Tasktr,Ytr,fun, params);

repelem(result_singleTask{kk},5)
kk = 9;
%% YMatirx should be a numData * numLines for direct accuracy
YMatrix1 =[result_simpleSum{kk};  result_covariance{kk}; result_CGPUCB{kk};repelem(result_singleTask{kk},5)]';

kk = 8;
%% YMatirx should be a numData * numLines for inference accuracy
YPre =[result_simpleSum{kk};  result_covariance{kk}; result_CGPUCB{kk};repelem(result_singleTask{kk},5)]';
YMatrix2 = zeros(size(YPre));
for jj=1:100
YMatrix2(jj,:) = max(YPre(1:jj,:),[],1);
end

% repelem(result_singleTask{kk},5) for each point, we evalute 5 times so
% the loss should also repeat 5 times
plot(result_simpleSum{kk})
hold on;
plot(result_covariance{kk})
hold on;
plot(result_CGPUCB{kk})
hold on;
plot(repelem(result_singleTask{kk},5))
hold off;


%% for repeated  experiments to get the average peroformance
% 10 times for first try
bounds = repmat([-3 3],2,1);

fun = @(x,i) Evaluate_Arcene(x, i, Data, Label,indices);
% fun([0 2],[ 1 2])
n =20; % should be even number
rng(0,'twister'); % For reproducibility
XtrRaw = boGetInitPts(bounds,n, []);
Tasktr = [repmat(1, n/5, 1); repmat(2, n/5, 1);repmat(3, n/5, 1);repmat(4, n/5, 1);repmat(5, n/5, 1)];
for i = 1:n
Ytr(i,1) = fun(XtrRaw(i,:),Tasktr(i));
end

tic
% main iterations
for iter = 1:20
    
params.numTasks = NumFold;
params.fmincon.gridSize = 10000;
params.bounds = repmat([-3 3],2,1);
params.BetaUpdate= 'simple';
% params.BetaUpdate= 'complex';
params.prob= 0.2; 
params.NumSearch = 1;
% params.Task2Numeric = [];
params.KernelFun = 'ardmatern52';
params.T = 100;
params.BatchSize = 1;
params.Task2Numeric = [0 1 1 1 1;1 0 1 1 1;1 1 0 1 1; 1 1 1 0 1; 1 1 1 1 0];
% params.Task2Numeric = [1 0 0;0 1 0;0 0 1];
params.Strategy = 'MT-GPUCB';
% params.Strategy = 'Random';
% params.optimalSolution = [0 0];
params.NewAcq = 'SimpleSum';
params.MultiEvaluation = 'On';


params.Seed = iter;


summary_simpleSum{iter} = myBO(XtrRaw,Tasktr,Ytr,fun, params);


clear params
disp(iter)
end

toc



%% prepare for simple regret data
for iter = 1:20
record_simpleSum(iter,:) = summary_simpleSum{iter}{9}; 
end

mean(record_simpleSum)
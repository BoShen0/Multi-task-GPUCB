%% Multiple tasks
params.numTasks = 3;
params.fmincon.gridSize = 10000;
params.bounds = repmat([-4 5],4,1);
params.BetaUpdate= 'simple';
% params.BetaUpdate= 'complex';
params.prob= 0.01; 
params.NumSearch = 1;
% params.Task2Numeric = [];
params.KernelFun = 'ardmatern52';
params.T = 60;
params.BatchSize = 1;
params.Task2Numeric = [10 5 2 10;
                       10 5 3 9;
                       8 6 2 10];
% params.Task2Numeric = [1 0 0;0 1 0;0 0 1];
params.Strategy = 'MT-GPUCB';
% params.Strategy = 'Random';
params.optimalSolution = [0 0 0 0];
params.NewAcq = 'SimpleSum';

fun = @(x,i) Evaluate_powell(x, i);

n =78; % should be even number
rng(0,'twister'); % For reproducibility
XtrRaw = boGetInitPts(params.bounds,n, []);
Tasktr = [repmat(1, n/3, 1); repmat(2, n/3, 1);
          repmat(3, n/3, 1)];
for i = 1:n
Ytr(i,1) = fun(XtrRaw(i,:),Tasktr(i));
end

% for simple sum 
result_simpleSum = myBO(XtrRaw,Tasktr,Ytr,fun, params);

% for covariance 
params.prob= 0.05; 
params.BetaUpdate= 'two';
params.NewAcq = 'Covariance';
result_covariance = myBO(XtrRaw,Tasktr,Ytr,fun, params);

% GPUCB
params.Strategy = 'CGP-UCB';
params.NewAcq = 'CGP-UCB';
result_CGPUCB = myBO(XtrRaw,Tasktr,Ytr,fun, params);

% single task
params.Strategy = 'MT-GPUCB';
params.NewAcq = 'Covariance';
params.BatchSize = 3;
params.T = 20;
result_singleTask = myBO(XtrRaw,Tasktr,Ytr,fun, params);

kk = 9;
%% YMatirx should be a numData * numLines for direct regret
YMatrix1 =[result_simpleSum{kk};  result_covariance{kk}; result_CGPUCB{kk};...
    repelem(result_singleTask{kk},3)]';

kk = 10;
%% YMatirx should be a numData * numLines for cumulative regret
Direct_regret = repelem(result_singleTask{7},3);
for i =1:60
    YY(i)=sum(Direct_regret(1:i));
end
YMatrix2 =[result_simpleSum{kk};  result_covariance{kk}; result_CGPUCB{kk};YY]';

kk = 8;
%% YMatirx should be a numData * numLines for inference accuracy
YPre =[result_simpleSum{kk};  result_covariance{kk}; result_CGPUCB{kk};repelem(result_singleTask{kk},3)]';
YMatrix3 = zeros(size(YPre));
for jj=1:60
YMatrix3(jj,:) = min(YPre(1:jj,:),[],1);
end


%% for repeated  experiments to get the average peroformance
% 10 times for first try
bounds = repmat([-4 5],4,1);

fun = @(x,i) Evaluate_powell(x, i);

n =78; % should be even number
rng(0,'twister'); % For reproducibility
XtrRaw = boGetInitPts(bounds,n, []);
Tasktr = [repmat(1, n/3, 1); repmat(2, n/3, 1);
          repmat(3, n/3, 1)];
for i = 1:n
Ytr(i,1) = fun(XtrRaw(i,:),Tasktr(i));
end


tic
% main iterations
for iter = 1:50
    
params.numTasks = 3;
params.fmincon.gridSize = 10000;
params.bounds = repmat([-4 5],4,1);
params.BetaUpdate= 'simple';
% params.BetaUpdate= 'complex';
params.prob= 0.01; 
params.NumSearch = 1;
% params.Task2Numeric = [];
params.KernelFun = 'ardmatern52';
params.T = 90;
params.BatchSize = 1;
params.Task2Numeric = [10 5 2 10;
                       10 5 3 9;
                       8 6 2 10];
% params.Task2Numeric = [1 0 0;0 1 0;0 0 1];
params.Strategy = 'MT-GPUCB';
% params.Strategy = 'Random';
params.optimalSolution = [0 0 0 0];
params.NewAcq = 'SimpleSum';

params.Seed = iter;


summary_simpleSum{iter} = myBO(XtrRaw,Tasktr,Ytr,fun, params);


clear params
disp(iter)
end

toc


%% prepare for simple regret data
for iter = 1:50
record_simpleSum(iter,:) = summary_simpleSum{iter}{9}; 
end

mean(record_simpleSum)

std(record_simpleSum)
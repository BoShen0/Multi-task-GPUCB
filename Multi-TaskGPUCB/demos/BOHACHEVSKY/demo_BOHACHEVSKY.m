%% for repeated  experiments to get the average peroformance
% 10 times for first try
bounds = repmat([-50 50],2,1);

fun = @(x,i) Evaluate_boha(x, i);

n =15; % should be even number
rng(0,'twister'); % For reproducibility
XtrRaw = boGetInitPts(bounds,n, []);
Tasktr = [repmat(1, n/3, 1); repmat(2, n/3, 1);repmat(3, n/3, 1)];
for i = 1:n
Ytr(i,1) = fun(XtrRaw(i,:),Tasktr(i));
end

tic
% main iterations
for iter = 1:50
    
params.numTasks = 3;
params.fmincon.gridSize = 10000;
params.bounds = repmat([-50 50],2,1);
params.BetaUpdate= 'simple';
% params.BetaUpdate= 'complex';
params.prob= 0.1; 
params.NumSearch = 1;
% params.Task2Numeric = [];
params.KernelFun = 'ardmatern52';
params.T = 30;
params.BatchSize = 1;
params.Task2Numeric = [1 2;0.8 1.8;1.2 2.2];
% params.Task2Numeric = [1 0 0;0 1 0;0 0 1];
params.Strategy = 'MT-GPUCB';
% params.Strategy = 'Random';
params.optimalSolution = [0 0];
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
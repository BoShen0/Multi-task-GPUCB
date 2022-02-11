params.numTasks = 3;
params.fmincon.gridSize = 10000;
params.bounds = repmat([-10 10],2,1);
params.BetaUpdate= 'simple';
% params.BetaUpdate= 'complex';
params.prob= 0.02; 
params.NumSearch = 1;
% params.Task2Numeric = [];
params.KernelFun = 'ardmatern52';
params.T = 20;
params.BatchSize = 1;
params.Task2Numeric = [0.9 1 1;1 0.9 1;1 1 0.9];
% params.Task2Numeric = [1 0 0;0 1 0;0 0 1];
params.Strategy = 'MT-GPUCB';
% params.Strategy = 'Random';
params.optimalSolution = [0 0];
params.NewAcq = 'Covariance';

fun = @(x,i) Evaluate_boha(x, i);

rng(0,'twister'); % For reproducibility
n =3;
XtrRaw = randi([-10 10], n,2);
XtrRaw = repmat(XtrRaw, 3,1);
for i=1:n
Ytr(i,1) = -boha1(XtrRaw(i,:));
Ytr(n+i,1) = -boha2(XtrRaw(i,:));
Ytr(2*n+i,1) = -boha3(XtrRaw(i,:));
end
% Ytr = -Ytr/30000;
Tasktr = [repmat(1, n, 1); repmat(2, n, 1); repmat(3, n, 1)];

[acq, uncerts, covariance] = acqMTGPUCB([2 2], GPmodel, 1, numTasks, Task2Numeric,XtrTask, 'Covariance', prob, BetaUpdate, params.bounds);
sqrt(diag(covariance)) - uncerts
mean(covariance,'all')
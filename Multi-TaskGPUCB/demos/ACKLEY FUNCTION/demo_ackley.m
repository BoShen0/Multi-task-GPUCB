[X,Y] = meshgrid(-32.768:0.6:32.768);                                
Z1= zeros(size(X));
for i=1:size(X,1)
    for j=1:size(X,2)
    Z1(i,j) = ackley([X(i,j),Y(i,j)]);
    end
end

Z2= zeros(size(X));
for i=1:size(X,1)
    for j=1:size(X,2)
    Z2(i,j) = ackley([X(i,j),Y(i,j)],15,0.3);
    end
end

   
figure, surf(X,Y,Z1)
title('Ackley Function 1');
view([-45.54 9.29964912280702]);
figure, surf(X,Y,Z2)
title('Ackley Function 2');
view([-45.54 9.29964912280702]);

%% Single task
% params.numTasks = 1;
% params.fmincon.gridSize = 10000;
% params.bounds = repmat([-32.768 32.768],2,1);
% params.BetaUpdate= 'simple';
% params.prob= 0.2; 
% params.NumSearch = 1;
% params.Task2Numeric = [];
% params.KernelFun = 'ardmatern52';
% params.T = 100;
% params.BatchSize = 1;
% % params.NumSearch = 10;
% params.optimalSolution = [0 0];
% 
% 
% 
% fun = @(x) -ackley(x); 
% n =10;
% rng(0,'twister'); % For reproducibility
% Xtr = boGetInitPts(params.bounds,n, []);
% for i = 1:n
% Ytr(i,1) = fun(Xtr(i,:));
% end
% 
% rere1 = myBO(Xtr,[],Ytr,fun, params); % KMPP inital
% rere2 = myBO(Xtr,[], Ytr,fun, params); % LHS initial
% kk = 10;
% plot(rere1{kk})
% hold on;
% plot(rere2{kk})
% hold off;
% 
%% for two tasks BO with same optimal solution
% This is the basic set up
% params.numTasks = 2;
% params.fmincon.gridSize = 10000;
% params.bounds = repmat([-32.768 32.768],2,1);
% params.BetaUpdate= 'simple';
% % params.BetaUpdate= 'complex';
% params.prob= 0.2; 
% params.NumSearch = 1;
% % params.Task2Numeric = [];
% params.KernelFun = 'ardmatern52';
% params.T = 30;
% params.BatchSize = 1;
% params.Task2Numeric = [0.2 0.2;0.15 0.3];
% % params.Task2Numeric = [1 0 0;0 1 0;0 0 1];
% params.Strategy = 'MT-GPUCB';
% % params.Strategy = 'Random';
% params.optimalSolution = [0 0];
% params.NewAcq = 'SimpleSum';
% n =10; % should be even number

params.numTasks = 2;
params.fmincon.gridSize = 10000;
params.bounds = repmat([-32.768 32.768],2,1);
params.BetaUpdate= 'simple';
% params.BetaUpdate= 'complex';
params.prob= 0.2; 
params.NumSearch = 1;
% params.Task2Numeric = [];
params.KernelFun = 'ardmatern52';
params.T = 30;
params.BatchSize = 1;
params.Task2Numeric = [0.2 0.2;0.15 0.3];
% params.Task2Numeric = [1 0 0;0 1 0;0 0 1];
params.Strategy = 'MT-GPUCB';
% params.Strategy = 'Random';
params.optimalSolution = [0 0];
params.Seed = 5;
params.NewAcq = 'SimpleSum';

fun = @(x,i) Evaluate_ackley(x, i);

n =10; % should be even number
rng(0,'twister'); % For reproducibility
XtrRaw = boGetInitPts(params.bounds,n, []);
Tasktr = [repmat(1, n/2, 1); repmat(2, n/2, 1)];
for i = 1:n
Ytr(i,1) = fun(XtrRaw(i,:),Tasktr(i));
end

% result = myBO(XtrRaw,Tasktr,Ytr,fun, params);



% params.prob= 0.02; 
% result_simpleSum2 = myBO(XtrRaw,Tasktr,Ytr,fun, params);
% 
% params.KernelFun = 'ardsquaredexponential';
% result_simpleSumArdSquareExp = myBO(XtrRaw,Tasktr,Ytr,fun, params);
% 
% params.KernelFun = 'ardrationalquadratic';
% result_simpleSumARDRationQuad = myBO(XtrRaw,Tasktr,Ytr,fun, params);


% for simple sum 
params.NewAcq = 'SimpleSum';
result_simpleSum = myBO(XtrRaw,Tasktr,Ytr,fun, params);

% for covariance 
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
params.BatchSize = 2;
params.T = 15;
result_singleTask = myBO(XtrRaw,Tasktr,Ytr,fun, params);

% fixed task
% params.Strategy = 'FixTask';
% params.NewAcq = 'Covariance';
% params.BatchSize = 1;
% params.T = 30;
% result_fixTask = myBO(XtrRaw,Tasktr,Ytr,fun, params);

%% for repeated  experiments to get the average peroformance
% 10 times for first try
bounds = repmat([-32.768 32.768],2,1);

fun = @(x,i) Evaluate_ackley(x, i);

n =10; % should be even number
rng(0,'twister'); % For reproducibility
XtrRaw = boGetInitPts(bounds,n, []);
Tasktr = [repmat(1, n/2, 1); repmat(2, n/2, 1)];
for i = 1:n
Ytr(i,1) = fun(XtrRaw(i,:),Tasktr(i));
end

tic
% main iterations
for iter = 1:50
params.numTasks = 2;
params.fmincon.gridSize = 10000;
params.bounds = repmat([-32.768 32.768],2,1);
params.BetaUpdate= 'simple';
% params.BetaUpdate= 'complex';
params.prob= 0.2; 
params.NumSearch = 1;
% params.Task2Numeric = [];
params.KernelFun = 'ardmatern52';
params.T = 60;
params.BatchSize = 1;
params.Task2Numeric = [0.2 0.2;0.15 0.3];
% params.Task2Numeric = [1 0 0;0 1 0;0 0 1];
params.Strategy = 'MT-GPUCB';
% params.Strategy = 'Random';
params.optimalSolution = [0 0];
params.Seed = iter; % rng(params.Seed+5*t)
params.NewAcq = 'SimpleSum';


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
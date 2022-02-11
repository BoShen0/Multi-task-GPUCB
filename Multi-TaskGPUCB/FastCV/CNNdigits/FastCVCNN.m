% This experiment is for fast validation method
% using CNN for digit classification
%% prepare data for fast cross validation
NumFold = 5;
rng(10)
indices = crossvalind('Kfold',2000,NumFold);

[XTrain,YTrain] = digitTrain4DArrayData;
X = XTrain(:,:,:,1:2000);
Y = YTrain(1:2000);


layers = [
    imageInputLayer([28 28 1])
    
    convolution2dLayer(3,8,'Padding','same')
    batchNormalizationLayer
    reluLayer   
    
    maxPooling2dLayer(2,'Stride',2)
    
    convolution2dLayer(3,16,'Padding','same')
    batchNormalizationLayer
    reluLayer   
    
    maxPooling2dLayer(2,'Stride',2)
    
    convolution2dLayer(3,32,'Padding','same')
    batchNormalizationLayer
    reluLayer   
    
    fullyConnectedLayer(10)
    softmaxLayer
    classificationLayer];
tic
[Y_eval]= Evaluate_CNNDigit([-1 -1 0.1], 1:5, X,Y,indices,layers);
toc

%% for two tasks BO with same optimal solution
% 10^nextPt(1) -- 'InitialLearnRate';10^nextPt(2) -- 'L2Regularization'; 
% nextPt(3) -- 'Momentum' 
% params.numTasks = NumFold;
% params.fmincon.gridSize = 10000;
% params.bounds = [-5 0;-5 0;0.4 1]; 
% params.BetaUpdate= 'simple';
% % params.BetaUpdate= 'complex';
% params.prob= 0.2; 
% params.NumSearch = 1;
% % params.Task2Numeric = [];
% params.KernelFun = 'ardmatern52';
% params.T = 50;
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
params.bounds = [-5 0;-5 0;0.4 1]; 
params.BetaUpdate= 'simple';
% params.BetaUpdate= 'complex';
params.prob= 0.2; 
params.NumSearch = 1;
% params.Task2Numeric = [];
params.KernelFun = 'ardmatern52';
params.T = 50;
params.BatchSize = 1;
params.Task2Numeric = [0 1 1 1 1;1 0 1 1 1;1 1 0 1 1; 1 1 1 0 1; 1 1 1 1 0];
% params.Task2Numeric = [1 0 0;0 1 0;0 0 1];
params.Strategy = 'MT-GPUCB';
% params.Strategy = 'Random';
% params.optimalSolution = [0 0];
params.NewAcq = 'SimpleSum';
params.MultiEvaluation = 'On';

fun = @(x,i) Evaluate_CNNDigit(x, i, X,Y,indices,layers);
% fun([0 2],[ 1 2])
n =20; % should be even number
rng(0,'twister'); % For reproducibility
XtrRaw = boGetInitPts(params.bounds,n, []);
Tasktr = [repmat(1, n/5, 1); repmat(2, n/5, 1);repmat(3, n/5, 1);repmat(4, n/5, 1);repmat(5, n/5, 1)];

tic
for i = 1:n
Ytr(i,1) = fun(XtrRaw(i,:),Tasktr(i));
disp(i)
end
toc

% for simple sum 
result_simpleSum = myBO(XtrRaw,Tasktr,Ytr,fun, params);

% for covariance
params.BetaUpdate= 'complex'; % this is used to show the performance
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
params.T = 10;
result_singleTask = myBO(XtrRaw,Tasktr,Ytr,fun, params);

repelem(result_singleTask{kk},5)
kk = 9;
%% YMatirx should be a numData * numLines for direct accuracy
YMatrix1 =[result_simpleSum{kk};  result_covariance{kk}; result_CGPUCB{kk};repelem(result_singleTask{kk},5)]';

kk = 8;
%% YMatirx should be a numData * numLines for inference accuracy
YPre =[result_simpleSum{kk};  result_covariance{kk}; result_CGPUCB{kk};repelem(result_singleTask{kk},5)]';
YMatrix2 = zeros(size(YPre));
for jj=1:50
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
bounds = [-5 0;-5 0;0.4 1];

fun = @(x,i) Evaluate_CNNDigit(x, i, X,Y,indices,layers);
% fun([0 2],[ 1 2])
n =20; % should be even number
rng(0,'twister'); % For reproducibility
XtrRaw = boGetInitPts(bounds,n, []);
Tasktr = [repmat(1, n/5, 1); repmat(2, n/5, 1);repmat(3, n/5, 1);repmat(4, n/5, 1);repmat(5, n/5, 1)];

tic
for i = 1:n
Ytr(i,1) = fun(XtrRaw(i,:),Tasktr(i));
disp(i)
end
toc


tic
% main iterations
for iter = 1:20
    
params.numTasks = NumFold;
params.fmincon.gridSize = 10000;
params.bounds = [-5 0;-5 0;0.4 1]; 
params.BetaUpdate= 'simple';
% params.BetaUpdate= 'complex';
params.prob= 0.2; 
params.NumSearch = 1;
% params.Task2Numeric = [];
params.KernelFun = 'ardmatern52';
params.T = 50;
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

for iter = 1:20
record_simpleSum(iter,:) = summary_simpleSum{iter}{9}; 
end

mean(record_simpleSum)
[X,Y] = meshgrid(-10:0.2:10);                                
Z1= zeros(size(X));
for i=1:size(X,1)
    for j=1:size(X,2)
    Z1(i,j) = levy([X(i,j),Y(i,j)],[1 10 1]);
    end
end

Z2= zeros(size(X));
for i=1:size(X,1)
    for j=1:size(X,2)
    Z2(i,j) = levy([X(i,j),Y(i,j)],[1.2 9 1.2]);
    end
end

   
figure, surf(X,Y,Z1)
title('Levy Function 1');

figure, surf(X,Y,Z2)
title('Levy Function 2');




%% Multiple tasks
params.numTasks = 2;
params.fmincon.gridSize = 10000;
params.bounds = repmat([-10 10],2,1);
params.BetaUpdate= 'simple';
% params.BetaUpdate= 'complex';
params.prob= 0.1; 
params.NumSearch = 1;
% params.Task2Numeric = [];
params.KernelFun = 'ardmatern52';
params.T = 40;
params.BatchSize = 1;
params.Task2Numeric = [1 1 1;1.2 0.9 1.2];
% params.Task2Numeric = [1 0 0;0 1 0;0 0 1];
params.Strategy = 'MT-GPUCB';
% params.Strategy = 'Random';
params.optimalSolution = [1 1];
params.NewAcq = 'SimpleSum';

fun = @(x,i) Evaluate_levy(x, i);

n =30; % should be even number
rng(0,'twister'); % For reproducibility
XtrRaw = boGetInitPts(params.bounds,n, []);
Tasktr = [repmat(1, n/2, 1); repmat(2, n/2, 1)];
for i = 1:n
Ytr(i,1) = fun(XtrRaw(i,:),Tasktr(i));
end

% for simple sum 
result_simpleSum = myBO(XtrRaw,Tasktr,Ytr,fun, params);

% for covariance 
params.BetaUpdate= 'complex';
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
params.T = 20;
result_singleTask = myBO(XtrRaw,Tasktr,Ytr,fun, params);

kk = 9;
%% YMatirx should be a numData * numLines for direct regret
YMatrix1 =[result_simpleSum{kk};  result_covariance{kk}; result_CGPUCB{kk};...
    repelem(result_singleTask{kk},2)]';

kk = 10;
%% YMatirx should be a numData * numLines for cumulative regret
Direct_regret = repelem(result_singleTask{7},2);
for i =1:40
    YY(i)=sum(Direct_regret(1:i));
end
YMatrix2 =[result_simpleSum{kk};  result_covariance{kk}; result_CGPUCB{kk};YY]';

kk = 8;
%% YMatirx should be a numData * numLines for inference accuracy
YPre =[result_simpleSum{kk};  result_covariance{kk}; result_CGPUCB{kk};repelem(result_singleTask{kk},2)]';
YMatrix3 = zeros(size(YPre));
for jj=1:40
YMatrix3(jj,:) = min(YPre(1:jj,:),[],1);
end
% ax = gca;
% outerpos = ax.OuterPosition;
% ti = ax.TightInset; 
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3);
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];
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
bounds = repmat([-10 10],2,1);

fun = @(x,i) Evaluate_levy(x, i);

n =30; % should be even number
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
params.bounds = repmat([-10 10],2,1);
params.BetaUpdate= 'simple';
% params.BetaUpdate= 'complex';
params.prob= 0.1; 
params.NumSearch = 1;
% params.Task2Numeric = [];
params.KernelFun = 'ardmatern52';
params.T = 40;
params.BatchSize = 1;
params.Task2Numeric = [1 1 1;1.2 0.9 1.2];
% params.Task2Numeric = [1 0 0;0 1 0;0 0 1];
params.Strategy = 'MT-GPUCB';
% params.Strategy = 'Random';
params.optimalSolution = [1 1];
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
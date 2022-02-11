% Copyright (c) 2020 Bo Shen
% This function is partially adapted from the codes for the following paper
% 1. Hernández-Lobato, José Miguel, Matthew W. Hoffman, and Zoubin Ghahramani. "Predictive entropy search for efficient global 
% optimization of black-box functions." Advances in neural information processing systems. 2014.
% 2. Wang, Zi, and Stefanie Jegelka. "Max-value entropy search for efficient Bayesian optimization." arXiv preprint arXiv:1703.01968 (2017).
% 3. Kandasamy, Kirthevasan, et al. "Gaussian process bandit optimisation with multi-fidelity evaluations." 
% Advances in Neural Information Processing Systems. 2016.
function results = myBO(Xtr, Tasktr, Ytr, fun, params)
% Input: Xtr - raw value for the initial training data; Tasktr - task index
% for Xtr; Ytr - response variable; fun - the evaluate function; params -
% all parameters used in our multi-task BO.
% fun: the function to be maximized
%% pre-processing for params
if ~isfield(params, 'bounds');  error('Unknown bounds'); end
if ~isfield(params, 'numTasks'); error('Unknown Number of Tasks'); end
if ~isfield(params, 'BetaUpdate'); params.BetaUpdate = 'simple'; end
if ~isfield(params, 'prob'); params.prob = 0.1; end
if ~isfield(params, 'NumSearch'); params.NumSearch = 1; end
if ~isfield(params, 'T'); params.T = 50; end
if ~isfield(params, 'KernelFun'); params.KernelFun = 'ardmatern52'; end
if ~isfield(params, 'KernelParameters_Initial'); params.KernelParameters_Initial = []; end
if ~isfield(params, 'Basis'); params.Basis = 'constant'; end
if ~isfield(params, 'BatchSize'); params.BatchSize = 1; end
if ~isfield(params, 'NewAcq'); params.NewAcq = 'Covariance'; end
if ~isfield(params, 'Seed'); params.Seed = 10; end

if ~isfield(params.fmincon, 'gridSize'); params.fmincon.gridSize = 10000; end

numTasks = params.numTasks;  % number of tasks
BetaUpdate = params.BetaUpdate;  % beta updating rule
prob = params.prob;   % prob for beta

T =  params.T; % number of iteration
KernelFun = params.KernelFun; % kernel function
% inital kernel parameters, which is required for customized kernel
% function
KernelParameters_Initial = params.KernelParameters_Initial;
Basis = params.Basis;  % GP mean function
BatchSize = params.BatchSize;  % number of task batch

if numTasks == 1
    params.Task2Numeric = [];
else
    if isempty(params.Task2Numeric)
        error('Undefined Task2Numeric');
    end
end
Task2Numeric = params.Task2Numeric; % task features

if numTasks>1 && ~isfield(params, 'Strategy')
    error('Unknown Task Strategy: it is required for numTasks>1');
end
%% set up required/necessary variables 
guesses = Xtr;
XtrTask = []; % default setting for numTasks == 1.
 
if isfield(params, 'optimalSolution')
        if numTasks == 1
           fOptimal = fun(params.optimalSolution);
        else
           fOptimal = fun(params.optimalSolution,1:numTasks);
        end
end
    
%% main iteration for MTBO
% rng(params.Seed)

for t=1:params.T
    
    tStart = tic; 

    [GPmodel] = MTGPRegression(Xtr, Tasktr, Ytr, ...
    Task2Numeric , KernelFun, KernelParameters_Initial, Basis);
%     KernelParameters_Initial = GPmodel.KernelInformation.KernelParameters;
    mufun = @(x) muEST(x, GPmodel, numTasks, Task2Numeric);
    [guessPt, ~] = globalMaximization(mufun, params.bounds(:,1), params.bounds(:,2),guesses,0);
    guesses = [guesses; guessPt];
    params.guesses = guesses;
    
    if numTasks > 1
    XtrTask = [Xtr Task2Numeric(Tasktr,:)];
    end
    
%   rng(10) this is the old one    
    rng(params.Seed+5*t) %Ackley
%     rng(params.Seed) % this is for the experiemnt purpose
    if strcmp(params.Strategy,'CGP-UCB')
        rng(params.Seed+5*t)
    end
    
    [nextPt, nextTasksBatch, ~] = strategyMTGPUCB_fmincon(t,...
    GPmodel, XtrTask, params);
    if numTasks == 1
        Xtr = [Xtr;nextPt];
        Ytr = [Ytr;fun(nextPt)];
    else
        Xtr = [Xtr;repmat(nextPt,BatchSize,1)];
        Tasktr = [Tasktr;nextTasksBatch];
        YYYVal = fun(nextPt,1:numTasks);
        Ytr = [Ytr;YYYVal(nextTasksBatch,:)];   
    end

% calculate the regret if exist params.optimal solution
    if isfield(params, 'optimalSolution')
        if numTasks == 1
            DirectRegret(t) = fOptimal - fun(nextPt);
            InferRegret(t) =  fOptimal - fun(guessPt);
        else
            DirectRegret(t) = sum(fOptimal -  YYYVal)/numTasks;
            InferRegret(t) = sum(fOptimal -fun(guessPt,1:numTasks))/numTasks;
        end
        SimpleRegret(t) = min(DirectRegret);
        SumRegret(t) = sum(DirectRegret);
    end
 
% finish calculation of the regret

% calculate the multi-task loss if  params.optimal solution does not exist
% and multi-evaluation exist
  if ~isfield(params, 'optimalSolution') && isfield(params, 'MultiEvaluation')
        DirectLoss(t) = sum(  YYYVal)/numTasks;
        InferLoss(t)  = sum(fun(guessPt,1:numTasks))/numTasks;
        MaxDirectLoss(t) = max(DirectLoss);
        MaxInferLoss(t) = max(InferLoss);
  end
% end loss calculation

    tEnd = toc(tStart);
    time_iter(t) =   tEnd;
    GuessPts(t,:) = guessPt;
    disp(['Iteration ' num2str(t)   ': tested ' num2str(nextPt) ...
       '; taskBatch ' num2str(nextTasksBatch') '; funval ' num2str(Ytr((end - BatchSize +1):end)') ...
       '; iter_time '    num2str(tEnd)])
   
end

results{1} = Xtr;
results{2} = Tasktr;
results{3} = Ytr;
results{4} = GuessPts;
results{5} = time_iter;
results{6} = GPmodel;

if isfield(params, 'optimalSolution')
results{7} = DirectRegret;
results{8} = InferRegret;
results{9} = SimpleRegret;
results{10} = SumRegret;
end

if ~isfield(params, 'optimalSolution') && isfield(params, 'MultiEvaluation')
results{7} = DirectLoss;
results{8} = InferLoss;
results{9} = MaxDirectLoss;
results{10} = MaxInferLoss;
end

end


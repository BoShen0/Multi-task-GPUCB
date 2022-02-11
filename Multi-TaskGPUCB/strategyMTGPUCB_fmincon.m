function [nextPt, nextTasksBatch, nextPtAcq] = strategyMTGPUCB_fmincon(t,...
    GPmodel, XtrTask, params)
% Copyright (c) 2020 Bo Shen
% this function calculate the next point and task to evaluate, namely,
% nextPt and nextTaskBatch, respectively. nextPtAcq: acquistion value for
% the nextPt
% Input: GPmodel - A trained Gaussian process regression model; params - a
% structure that contains all the related components; XtrTask - The task
% index; t -  the number of iteration
 %% set up parameters
  numDims = size(params.bounds, 1);
  bounds = params.bounds;
  gridSize = params.fmincon.gridSize;
  prob =  params.prob;
  BetaUpdate = params.BetaUpdate;
  NumSearch = params.NumSearch;
  numTasks = params.numTasks;
  NewAcq = params.NewAcq;
  % generate random samples for optimiation
    Xgrid = rand_sample_interval(bounds(:,1), bounds(:,2), gridSize);
    if isfield(params, 'guesses')
    Xgrid = [Xgrid;params.guesses];
    end
%% main algorithm for the calculation   
if params.numTasks == 1
    y = arrayfun(@(ROWIDX) acqMTGPUCB(Xgrid(ROWIDX,:),GPmodel, t, 1, [], XtrTask,NewAcq, prob, BetaUpdate, bounds),...
    (1:gridSize).');
    [ySort,xIndex] = sort(y,'descend');
    maxVal = ySort(1);
    % our problem is maximization problem, fmincon is for minization problem,
    % so negative our object function
    neg_target = @(x) -acqMTGPUCB(x,GPmodel, t, 1, [],XtrTask,NewAcq, prob, BetaUpdate, bounds);
    gradObj = 'off';
    startMultiple = Xgrid(xIndex(1:NumSearch),:);
    
    optimum = zeros(NumSearch, numDims);
    fval = zeros(NumSearch, 1);
    
    for i=1:NumSearch
        [xval, functionval] = fmincon(neg_target, startMultiple(i,:), [], [], [], [], bounds(:,1), bounds(:,2), [], ...
        optimset('MaxFunEvals', 1000, 'TolX', eps, 'Display', 'off', 'GradObj', gradObj));
        functionval = -functionval;
        optimum(i,:)=  xval;
        fval(i,1) = functionval;
        if functionval < maxVal
        % In some situations, the optimization can fail. For example, when the
        % guesses contains the global optimum.
        %disp('fmincon in globalMaximization failed to return a value better than the initialization.')
            [xval2, functionval2] = fminsearch(neg_target, startMultiple(i,:));
            functionval2 = -functionval2;
            if functionval2 < maxVal
                %disp('fminsearch in globalMaximization failed to return a value better than the initialization.')
                fval(i,1) = maxVal;
                optimum(i,:) = startMultiple(1,:);
            else
                flag = 0;
                for ii = 1:numDims
                    if xval2(ii) < bounds(ii,1) || xval2(ii) >  bounds(ii,2)
                        flag = 1;
                        break;
                    end
                end
                if flag
                    %disp('fminsearch in globalMaximization failed to search within [xmin, xmax].')
                    fval(i,1) = maxVal;
                    optimum(i,:) = startMultiple(1,:);
                else
                    fval(i,1) = functionval2;
                    optimum(i,:) = xval2;                  
                end
            end  % functionval2 < maxVal
        end   % functionval < maxVal
    end   % NumSearch
    [nextPtAcq, maxIdx] = max(fval);
    nextPt = optimum(maxIdx,:);
    nextTasksBatch = 1;   
else   
    y = arrayfun(@(ROWIDX) acqMTGPUCB(Xgrid(ROWIDX,:),GPmodel, t, numTasks, params.Task2Numeric,XtrTask,NewAcq, prob, BetaUpdate, bounds),...
    (1:gridSize).');
    [ySort,xIndex] = sort(y,'descend');
    maxVal = ySort(1);
    % our problem is maximization problem, fmincon is for minization problem,
    % so negative our object function
    neg_target = @(x) -acqMTGPUCB(x,GPmodel, t, numTasks, params.Task2Numeric,XtrTask,NewAcq, prob, BetaUpdate, bounds);
    gradObj = 'off';
    startMultiple = Xgrid(xIndex(1:NumSearch),:);
    fval = zeros(NumSearch, 1);
    optimum = zeros(NumSearch, numDims);

    for i=1:NumSearch
        [xval, functionval] = fmincon(neg_target, startMultiple(i,:), [], [], [], [], bounds(:,1), bounds(:,2), [], ...
        optimset('MaxFunEvals', 1000, 'TolX', eps, 'Display', 'off', 'GradObj', gradObj));
        functionval = -functionval;
        optimum(i,:)=  xval;
        fval(i,1) = functionval;
        if functionval < maxVal
        % In some situations, the optimization can fail. For example, when the
        % guesses contains the global optimum.
        %disp('fmincon in globalMaximization failed to return a value better than the initialization.')
            [xval2, functionval2] = fminsearch(neg_target, startMultiple(i,:));
            functionval2 = -functionval2;
            if functionval2 < maxVal
                %disp('fminsearch in globalMaximization failed to return a value better than the initialization.')
                fval(i,1) = maxVal;
                optimum(i,:) = startMultiple(1,:);
            else
                flag = 0;
                for ii = 1:numDims
                    if xval2(ii) < bounds(ii,1) || xval2(ii) >  bounds(ii,2)
                        flag = 1;
                        break;
                    end
                end
                if flag
                    %disp('fminsearch in globalMaximization failed to search within [xmin, xmax].')
                    fval(i,1) = maxVal;
                    optimum(i,:) = startMultiple(1,:);
                else
                    fval(i,1) = functionval2;
                    optimum(i,:) = xval2;                  
                end
            end  % functionval2 < maxVal
        end   % functionval < maxVal
    end   % NumSearch
    [nextPtAcq, maxIdx] = max(fval);
    nextPt = optimum(maxIdx,:);     
    
    switch params.Strategy
    case 'MT-GPUCB'
    % First maximise the MF-GP-UCB acquisition function.
    [~, uncerts, ~] = acqMTGPUCB(nextPt,GPmodel, t, numTasks, params.Task2Numeric,XtrTask,NewAcq, prob, BetaUpdate, bounds);
    % Now determine the task batch
    [~, TaskRanks] = sort(uncerts,'descend');
    nextTasksBatch = TaskRanks(1:params.BatchSize);
    
    case 'CGP-UCB'
    nextTasksBatch = mod(t,numTasks)+1;
    
    case 'Random'
        rng('shuffle')
    % Now determine the task batch
    nextTasksBatch = randperm(numTasks,params.BatchSize)';

    case 'FixTask'    
    % Now determine the task batch
    nextTasksBatch = (1:params.BatchSize)';

    otherwise
    error('Unknown Strategy');
    end
end  % params.numTasks == 1

end


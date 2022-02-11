function [nextPt, nextTasksBatch, nextPtAcq] = strategyMTGPUCB(t,...
    funcHs, bounds, numTasks, prob, params)

  % Preliminaries and preprocessing
%   numFidels = numel(funcHs);
  numDims = size(bounds, 1);

   switch params.Strategy
  
       case 'MT-GPUCB'
  % First maximise the MF-GP-UCB acquisition function.
  acquisition = @(arg) acqMTGPUCB(arg, funcHs, t, numTasks, prob);
  [nextPtAcq, nextPt] = ...
    diRectWrap( acquisition, bounds, params.diRectParams);
  [nextPtAcq, uncerts] = acquisition(nextPt);

  % Now determine the task batch
[~, TaskRanks] = sort(uncerts,'descend');
nextTasksBatch = TaskRanks(1:params.BatchSize);

      case 'Random'
  acquisition = @(arg) acqMTGPUCB(arg, funcHs, t, numTasks, prob);
  [nextPtAcq, nextPt] = ...
    diRectWrap( acquisition, bounds, params.diRectParams);
  [nextPtAcq, uncerts] = acquisition(nextPt);
  % Now determine the task batch
  nextTasksBatch = randperm(numTasks,params.BatchSize);
  
      case 'FixTask'    
  acquisition = @(arg) acqMTGPUCB(arg, funcHs, t, params.BatchSize, prob);
  [nextPtAcq, nextPt] = ...
    diRectWrap( acquisition, bounds, params.diRectParams);
  [nextPtAcq, uncerts] = acquisition(nextPt);
  % Now determine the task batch
  nextTasksBatch = 1:params.BatchSize;
  
      otherwise
      error('Unknown Strategy');
   end
  

end


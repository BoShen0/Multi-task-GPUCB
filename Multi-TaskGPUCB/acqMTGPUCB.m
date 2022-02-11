function [acq, uncerts, covariance] = acqMTGPUCB(x, GPmodel, t, numTasks, Task2Numeric,XtrTask, NewAcq, prob, BetaUpdate, bounds)
% Copyright (c) 2020 Bo Shen
% This function calculates the acquistion function for point x under
% Gaussian process model, i.e., GPmodel, at iteration t. 
% x should be a row vector, current version can not support matrix input x
% Input: we need number of tasks (numTasks); task encoding matrix
% (Task2Numeric), which can be [] for numTasks==1; XtrTask, the traning
% data with task inforatmation, that will be used when NewAcq ==
% 'Covariance'; prob, BetaUpdate, and bounds are used to determine the
% Beta, which is the tradeoff term in GP-UCB type of algorithms.
% Output: acquistion value (acq); standard deviation (uncerts); covariance
% matirx of x @ different tasks, [] when numTasks == 1.

  numDims = size(x, 2);  % the dimension of the input
%   numTasks = numel(funcHs);
%   if ~isa(funcHs, 'cell') 
%     funcHs = {funcHs};
%   end'
%% determine the selection of beta 
   switch BetaUpdate  
       case 'two'
      beta_t = 2;
       case 'simple'
      beta_t = prob * numDims * log(2*numDims*t);
       case 'complex' 
      beta_t = 2*log(t^2*2*pi^2/(3*prob)) + 2*numDims*log(t^2*...
            numDims*max(bounds(:,2)-bounds(:,1))*(log(4*numDims/prob))^0.5);
       otherwise
      error('Donot know how to update beta');
   end
%   beta_t = 0.2 * numDims * log(2*numDims*numFidels*t); need to check the
%   definition
%   
  
%   beta_t = 2;
%% main iteration for calculation for the acquistion function
covariance = [];
  Sigma0 = GPmodel.Sigma;
if numTasks ==1
    [mu, sigma,~] = predict(GPmodel,x);
    uncerts = sqrt(sigma^2- Sigma0^2); % remove the observation noise variance
    acq = mu + sqrt(beta_t) * uncerts;
else
  indUCBs = zeros(numTasks, 1);
  uncerts = zeros(numTasks, 1);
  Mus = zeros(numTasks, 1);
  for i = 1:numTasks
 % get the mu and sigma for each task with the same x
    [mu, sigma,~] = predict(GPmodel,[x Task2Numeric(i,:)]); 
    Mus(i) = mu;
    uncerts(i) =  sqrt(sigma^2- Sigma0^2); % remove the observation noise variance
    indUCBs(i) = mu + sqrt(beta_t) * uncerts(i);
  end
  acq = sum(indUCBs);
  %% new acquistion function
  if strcmp(NewAcq,'Covariance')
      % call the model and calculate covariance matrix
      kfcn = GPmodel.Impl.Kernel.makeKernelAsFunctionOfXNXM(GPmodel.Impl.ThetaHat);
      Xte = [repmat(x,numTasks, 1) Task2Numeric];
      Ktrtr = kfcn(XtrTask, XtrTask) + diag(GPmodel.Sigma^2 * ones(size(XtrTask,1), 1));
      L = stableCholesky(Ktrtr);
      covariance = kfcn(Xte,Xte) - kfcn(Xte,XtrTask) * (L' \ (L \ kfcn(XtrTask,Xte)));
      acq = mean(Mus) + sqrt(beta_t) * sqrt(mean(covariance,'all'));
  elseif strcmp(NewAcq,'CGP-UCB')
  acq = indUCBs(mod(t,numTasks) + 1);
  end
end

end


function [GPmodel] = MTGPRegression(XtrRaw, Tasktr, Ytr, ...
    Task2Numeric, KernelFun, KernelParameters_Initial, Basis)
% Copyright (c) 2020 Bo Shen
% this function output the multi-task guassian process model
% Input: XtrRaw - training data without task information; Tasktr - task index
% for each XtrRaw; Ytr - response value; Task2Numeric - task encoding
% matrix; KernelFun - kernel function for fitting; KernelParameters_Initial
% - initial value for estmating the kernel parameters, which is required
% for customized kernel function; Basis - the basis of mean function in 
% gaussian process regression


% get the training data with numerical task information
numTrData = size(XtrRaw, 1);

if isempty(Task2Numeric)
      Xtr = XtrRaw;
else
    Numeric_tr = zeros(numTrData, size(Task2Numeric, 2));
    for ii = 1:numTrData
        Numeric_tr(ii,:) = Task2Numeric(Tasktr(ii),:);
    end    
    Xtr = [XtrRaw Numeric_tr];
end


% get the test data with numerical task information
  if numTrData == 0
    Ytr = zeros(0, 1);
  end
  
%% learn the GP
% mean functions we can explore, we should also know what is the
% optimization algorithm for estimating the hyper-parameters of the
% kernel function
if isempty(KernelParameters_Initial)
      GPmodel = fitrgp(Xtr,Ytr,'Basis',Basis, 'FitMethod','exact',...
      'KernelFunction',KernelFun);
else
      GPmodel = fitrgp(Xtr,Ytr,'Basis',Basis, 'FitMethod','exact',...
      'KernelFunction',KernelFun,'KernelParameters',KernelParameters_Initial);
end
    

%   Learned_KernelParameters = GPmodel.KernelInformation.KernelParameters;
%   Learned_noiseVar = GPmodel.Sigma;
%   Made a mistask for not putting a ^2 for the standard deviation
%   Ktrtr = KernelFun(Xtr, Xtr, Learned_KernelParameters) + diag(Learned_noiseVar^2 * ones(numTrData, 1));
%   Y_ = Ytr - meanFunc(Tasktr);
%   L = stableCholesky(Ktrtr);
%   alpha = L' \ (L \ Y_);
  
 
  % obtain the function handle
%   funcH = @(X, Task) GPComputeOutputs(X,Task, XtrRaw, Tasktr, Task2Numeric,...
% L, alpha, KernelFun,  Learned_KernelParameters, meanFunc);
  % Compute outputs for the test data
%   if ~isempty(Xte)
%     [teMean, teK, teStd] = funcH(Xte, Taskte);
%   else
%     teMean = []; teK = []; teStd = [];
%   end

end





function [yMu, yK, yStd] = GPComputeOutputs(XteRaw,Taskte, XtrRaw, Tasktr,Task2Numeric,...
L, alpha, KernelFun,  Learned_KernelParameters, meanFunc)

% get the training data with numerical task information
numTrData = size(XtrRaw, 1);
Numeric_tr = zeros(numTrData, size(Task2Numeric, 2));
for ii = 1:numTrData
    Numeric_tr(ii,:) = Task2Numeric(Tasktr(ii),:);
end    
Xtr = [XtrRaw Numeric_tr];

% get the test data with numerical task information
if ~isempty(XteRaw)
    Numeric_te = zeros(size(XteRaw, 1), size(Task2Numeric, 2));
    for ii = 1:size(XteRaw, 1)
        Numeric_te(ii,:) = Task2Numeric(Taskte(ii),:);
    end    
    Xte = [XteRaw Numeric_te];
else
    Xte=[];
end

  meanXte = meanFunc(Taskte);
  Ktetr = KernelFun(Xte, Xtr, Learned_KernelParameters);
  Ktete = KernelFun(Xte, Xte, Learned_KernelParameters);

  % Predictive Mean
  yMu = meanXte + Ktetr * alpha;
  % Predictive Variance
  V = L \ (Ktetr)';
  yK = Ktete - V'*V;
  yStd = sqrt(real(diag(yK)));
  
end


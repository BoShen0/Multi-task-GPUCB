function output = paperTan(Xtr, Tasktr, Ytr, fun, Tanparams)
% source: Tan, Matthias HY. "Bayesian optimization of expected quadratic 
% loss for multiresponse computer experiments with internal noise." 
% SIAM/ASA Journal on Uncertainty Quantification 8.3 (2020): 891-925.
% !!!!!!!the problem is in the minimization form

% clear all
% set up the global varibles for the implementation
format long g
global D11 thetaopt11 alphahat11 sigmahat11 GammaDDI11 Beta11 Y11 
global D12 thetaopt12 alphahat12 sigmahat12 GammaDDI12 Beta12 Y12
global D13 thetaopt13 alphahat13 sigmahat13 GammaDDI13 Beta13 Y13
global D14 thetaopt14 alphahat14 sigmahat14 GammaDDI14 Beta14 Y14
global optx11 optx13  
global dx target target2 QuadPoints P W Ohm m q tau11 tau12 tau14
%checked$
%% experimental setup
searchSize=Tanparams.searchSize; %\rho
noiselevel=Tanparams.noiselevel; %\rho
nadd=Tanparams.T; 
nodesigns=Tanparams.Repetition;
Taskinfo = Tanparams.Task2Numeric;
QuadPoints = Tanparams.Task2Numeric; % task features should be transfromed in advance
m0=size(QuadPoints,1); m=m0;%m=m0^2;

targetset=Tanparams.targetset;
if isfield(Tanparams, 'optimalSolution')
   optimalx = Tanparams.optimalSolution;
   optimalf = fun(optimalx,1:m0);
end

OriginalBounds = Tanparams.bounds;
NewBounds = repmat([0 1],size(OriginalBounds,1),1);
No=Inf;
nostart=1; 
%%%%%%%%%% 10^? control the convergence rate
if dx<5
options=optimoptions(@patternsearch,'MaxIter',10000,'TolFun',10^-8);
else
options=optimoptions(@patternsearch,'MaxIter',10000,'TolFun',10^-7);
end

n0=size(Xtr,1);


dx=size(Xtr,2); dz=size(QuadPoints,2); q=2;
W=[1];% W=[1 0; 0 1];

P=diag(ones(1,size(QuadPoints,1))/(2*sqrt(size(QuadPoints,1))));
% P=diag([(18-sqrt(30))/36 (18+sqrt(30))/36 (18+sqrt(30))/36 (18-sqrt(30))/36]/2);
%P=diag([(322-13*sqrt(70))/900 (322+13*sqrt(70))/900 128/225 (322+13*sqrt(70))/900 (322-13*sqrt(70))/900]/2);
% P=kron(P,P);
Ohm=kron(P,W);

if ~isfield(Tanparams, 'gridSize'); Tanparams.gridSize = 80000; end

if (searchSize+1)^dx <Tanparams.gridSize
gridx0=(fullfact( (searchSize+1)*ones(1,dx))-1)/searchSize; ngridx0=size(gridx0,1);
% gridx2=(fullfact((searchSize+1)*ones(1,dx))-1)/searchSize; ngridx2=size(gridx2,1); 
gridx3=(fullfact( (searchSize+1)*ones(1,dx))-1)/searchSize; 
else
    rng(10,'twister'); % For reproducibility
gridx0=lhsdesign(Tanparams.gridSize,dx); ngridx0=size(gridx0,1);
gridx3=lhsdesign(Tanparams.gridSize,dx);
end    
    
% [0.0133379278857657 1.99432496185344;
%            13.7858976716857          16.4466194466848];

% for i=1:nodesigns
%     Designs{i}(:,(dx+1):(dx+dz))=TruncNormInv(Designs{i}(:,(dx+1):(dx+dz)),noiselevel);
%     Designs2{i}(:,(dx+1):(dx+dz))=TruncNormInv(Designs2{i}(:,(dx+1):(dx+dz)),noiselevel);  
% end
Nrun=size(targetset,2)*nodesigns; 
targetset=kron(targetset,ones(1,nodesigns));

target=targetset(:,1); target2=repmat(target,m,1);
% Dev =  fun(optimalx,1:m0) -repmat(target,1,m);
% optimalEQL = sum(sum(Dev.*(W*Dev*P)));

Performance11=cell(1,Nrun); Performance12=cell(1,Nrun); Performance13=cell(1,Nrun); Performance14=cell(1,Nrun); 
% Distance11=zeros(nadd,Nrun); Distance12=zeros(nadd,Nrun); Distance13=zeros(nadd,Nrun); Distance14=zeros(nadd,Nrun); 
% for i=1:Nrun
%     Performance11{i}=zeros(nadd,4);
%     Performance12{i}=zeros(nadd,4);
%     Performance13{i}=zeros(nadd,4);      
%     Performance14{i}=zeros(nadd,4);          
% %     Performance2{i}=zeros(1,4);    
% end
%% Iterations for multiple runs
for run=1:Nrun  
    rng(10*run);
    rho = unifrnd(noiselevel,2*noiselevel);
    QuadPoints =  TruncNormInv(Taskinfo,rho);
    gridz=QuadPoints; ngridz=size(gridz,1);
   
n=n0;   
XtrNew = inputTransformer(Xtr, OriginalBounds, NewBounds);
D11=[XtrNew,QuadPoints(Tasktr,:)];
D12=D11;
D13=D11;
D14=D11;

Y11=Ytr';
Y12=Y11;
Y13=Y11;
Y14=Y11;

stop=0; count=0; EI11=[]; EI12=[]; EI13=[]; EI14=[];
estoptx11=zeros(nadd,dx); estoptx12=zeros(nadd,dx); estoptx13=zeros(nadd,dx); estoptx14=zeros(nadd,dx);

[thetaopt11,alphahat11,sigmahat11,GammaDDI11,Beta11]=GPFitXZV3(D11,Y11);
thetaopt12=thetaopt11; thetaopt13=thetaopt11; thetaopt14=thetaopt11; 
alphahat12=alphahat11; alphahat13=alphahat11; alphahat14=alphahat11; 
sigmahat12=sigmahat11; sigmahat13=sigmahat11; sigmahat14=sigmahat11;
GammaDDI12=GammaDDI11; GammaDDI13=GammaDDI11; GammaDDI14=GammaDDI11;
Beta12=Beta11; Beta13=Beta11; Beta14=Beta11; 

Q11=compEQL(D11,D11(:,1:dx),QuadPoints,P,W,target,thetaopt11,alphahat11,sigmahat11,GammaDDI11,Beta11);
[minQ11,index]=min(Q11);
minQ12=minQ11; minQ14=minQ11; 
Q0=compEQL(D14,gridx3,QuadPoints,P,W,target,thetaopt14,alphahat14,sigmahat14,GammaDDI14,Beta14);
[minQ0,index0]=min(Q0);
if(minQ0<minQ14)
    start=gridx3(index0,:);
else
    start=D14(index,1:dx);
end
EQL14=@(x)compEQL(D14,x,QuadPoints,P,W,target,thetaopt14,alphahat14,sigmahat14,GammaDDI14,Beta14);        
[estoptx,estobjval14,exitflag]=patternsearch(EQL14,start,[],[],[],[],zeros(1,dx),ones(1,dx),[],options);         
%% Main iterations for Bayesian optimization
while(stop==0)
count=count+1;    
 tStart = tic; 
 
    tic
    tau11=minQ11; 
    gridx=[gridx0; D11(:,1:dx); [estoptx; estoptx11(1:(count-1),:)]]; 
    ComputeEI11=@(x)CompEI(x,D11,thetaopt11,alphahat11,sigmahat11,GammaDDI11,Beta11,target2,QuadPoints,P,W,m,tau11);   
    if(count==1)
        Baseline=-ComputeEI11(estoptx);
        Decision=CompareEQL(D11,gridx,QuadPoints,P,W,target,thetaopt11,alphahat11,sigmahat11,GammaDDI11,Beta11,tau11,Baseline);        
    else
        Baseline=-ComputeEI11(estoptx11(count-1,:));
        Decision=CompareEQL(D11,gridx,QuadPoints,P,W,target,thetaopt11,alphahat11,sigmahat11,GammaDDI11,Beta11,tau11,Baseline);
    end
    gridx=gridx(Decision,:); ngridx=size(gridx,1);
    storecriterion11=zeros(ngridx,1); 
    Time =toc;
disp(['Run ' num2str(run)   ': Iteration ' num2str(count) ...
       '--- Module1 sub1 ' num2str(Time) '--- gridsize ' num2str(ngridx)])
   
tic
parfor i=1:ngridx %par
    storecriterion11(i)=ComputeEI11(gridx(i,:));
end
Time =toc;
disp(['Run ' num2str(run)   ': Iteration ' num2str(count) ...
       '--- Module1 sub2      ' num2str(Time) '--- gridsize ' num2str(ngridx)])
   
tic
[~,i0]=sort(storecriterion11); 
optx=zeros(nostart,dx); fval=zeros(nostart,1);
% parfor j=1:nostart
j=1;
    [optx(j,:),fval(j),exitflag] = patternsearch(ComputeEI11,gridx(i0(j),:),[],[],[],[],zeros(1,dx),ones(1,dx),[],options);
% end
[mEI,JJ]=min(fval); optx11=optx(JJ,:);
EI11=[EI11; -mEI];
Time =toc;
disp(['Run ' num2str(run)   ': Iteration ' num2str(count) ...
       '--- Module2 ' num2str(Time)])

% if((count==No)&&(run==1))
%     storecriterion11=zeros(ngridx0,1); 
% parfor i=1:ngridx0
%     storecriterion11(i)=ComputeEI11(gridx0(i,:));
% end    
% figure(501),contour([0:0.05:1]',[0:0.05:1]',reshape(storecriterion11(1:ngridx0),21,21)','ShowText','on')
% hold on
% plot(optx11(1),optx11(2),'^')
% input('press enter')
% end
tic
[EVar,i0]=CompPV(D11,thetaopt11,alphahat11,sigmahat11,GammaDDI11,Beta11,optx11,target,QuadPoints,m,P,W);
optz11=QuadPoints(i0,:);
D11=[D11; optx11 optz11];
fValue = fun(inputTransformer(optx11,NewBounds,OriginalBounds), 1:m0);
Y11=[Y11 fValue(i0)];
if isfield(Tanparams, 'optimalSolution')
Data_save= [inputTransformer(optx11,NewBounds,  OriginalBounds)  i0 fValue(i0) mean(fValue-optimalf)];
else
Data_save= [inputTransformer(optx11,NewBounds,  OriginalBounds)  i0 fValue(i0) mean(fValue)];
end
n=n+1;

[thetaopt11,alphahat11,sigmahat11,GammaDDI11,Beta11]=GPFitXZV3(D11,Y11,thetaopt11,alphahat11,sigmahat11);
Time =toc;
disp(['Run ' num2str(run)   ': Iteration ' num2str(count) ...
       '--- Module3 ' num2str(Time)])

% );,thetaopt11,alphahat11,sigmahat11
% [thetaopt12,alphahat12,sigmahat12,GammaDDI12,Beta12]=GPFitXZV3(D12,Y12,thetaopt12,alphahat12,sigmahat12);
% [thetaopt13,alphahat13,sigmahat13,GammaDDI13,Beta13]=GPFitXZV3(D13,Y13,thetaopt13,alphahat13,sigmahat13);
% [thetaopt14,alphahat14,sigmahat14,GammaDDI14,Beta14]=GPFitXZV3(D14,Y14,thetaopt14,alphahat14,sigmahat14);
tic
EQL11=@(x)compEQL(D11,x,QuadPoints,P,W,target,thetaopt11,alphahat11,sigmahat11,GammaDDI11,Beta11);
Q11=EQL11(D11(:,1:dx));
[minQ11,index]=min(Q11);
candidates=[gridx3; estoptx11(1:(count-1),:)];
Q0=EQL11(candidates);
[minQ0,index0]=min(Q0);
if(minQ0<minQ11)
    start=candidates(index0,:);
else
    start=D11(index,1:dx);
end
[estoptx11(count,:),estobjval11,exitflag]=patternsearch(EQL11,start,[],[],[],[],zeros(1,dx),ones(1,dx),[],options); 
Time =toc;
disp(['Run ' num2str(run)   ': Iteration ' num2str(count) ...
       '--- Module4 ' num2str(Time) ])

% if((count==No)&&(run==1))
% Objval=zeros(121,1);
% gx=(fullfact([11 11])-1)/10;
% for i=1:121
%     Objval(i)=EQL11(gx(i,:));
% end
% figure(1001),contour([0:0.1:1]',[0:0.1:1]',reshape(Objval,11,11)','ShowText','on')
% hold on
% plot(estoptx11(count,1),estoptx11(count,2),'^')
% [minObjval,I]=min(Objval);
% [gx(I,:) minObjval; estoptx11(count,:) estobjval11];
% input('press enter')
% end
% trueobjval11=TrueEQL(estoptx11(count,:));
% estoptx11V2=D11(index,1:dx); estobjval11V2=minQ11; trueobjval11V2=TrueEQL(estoptx11V2);

% Performance11{run}(count,:)=[(estobjval11-trueobjval11) (trueobjval11-optimalEQL) (estobjval11V2-trueobjval11V2) (trueobjval11V2-optimalEQL)];
% Performance12{run}(count,:)=[(estobjval12-trueobjval12) (trueobjval12-optimalEQL) (estobjval12V2-trueobjval12V2) (trueobjval12V2-optimalEQL)];
% Performance13{run}(count,:)=[(estobjval13-trueobjval13) (trueobjval13-optimalEQL) (estobjval13V2-trueobjval13V2) (trueobjval13V2-optimalEQL)];
% Performance14{run}(count,:)=[(estobjval14-trueobjval14) (trueobjval14-optimalEQL) (estobjval14V2-trueobjval14V2) (trueobjval14V2-optimalEQL)];
tEnd = toc(tStart);
Performance11{run}(count,:)= [tEnd Data_save];

disp(['Run ' num2str(run)   ': Iteration ' num2str(count) ...
       '--- time ' num2str(tEnd) '--- gridsize ' num2str(ngridx)])
disp(' ');   
if(count==nadd)
    stop=1;
end
end

%% save data for each run
end
output = Performance11;
end
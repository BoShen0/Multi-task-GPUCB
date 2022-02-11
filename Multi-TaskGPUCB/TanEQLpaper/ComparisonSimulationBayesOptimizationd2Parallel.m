function ComparisonSimulationBayesOptimizationd2Parallel
%checked$27042019
clear all
format long g
global D11 thetaopt11 alphahat11 sigmahat11 GammaDDI11 Beta11 Y11 
global D12 thetaopt12 alphahat12 sigmahat12 GammaDDI12 Beta12 Y12
global D13 thetaopt13 alphahat13 sigmahat13 GammaDDI13 Beta13 Y13
global D14 thetaopt14 alphahat14 sigmahat14 GammaDDI14 Beta14 Y14
global optx11 optx13  
global dx target target2 QuadPoints P W Ohm m q tau11 tau12 tau14
%checked$
No=Inf; 
noiselevels=[0.1]; lnoiselevels=length(noiselevels); nostart=1; options=optimoptions(@patternsearch,'MaxIter',10000,'TolFun',10^-9);
for NL=1:lnoiselevels
n0=10; nfinal=30; dxz=4; nadd=nfinal-n0; nodesigns=20;
load('TwentyDesigns,n0=10,nfinal=30.mat')

dx=dxz/2; dz=dxz/2; q=2;
W=[1 0; 0 1];
noiselevel=noiselevels(NL);

QuadPoints=([-sqrt(3/7+2/7*sqrt(6/5)) -sqrt(3/7-2/7*sqrt(6/5)) sqrt(3/7-2/7*sqrt(6/5)) sqrt(3/7+2/7*sqrt(6/5))]+1)'/2;
%QuadPoints=([-sqrt(5+2*sqrt(10/7))/3 -sqrt(5-2*sqrt(10/7))/3 0 sqrt(5-2*sqrt(10/7))/3 sqrt(5+2*sqrt(10/7))/3]+1)'/2;
m0=size(QuadPoints,1); m=m0^2;
QuadPoints=TruncNormInv([repmat(QuadPoints,m0,1) kron(QuadPoints,ones(m0,1))],noiselevel);
P=diag([(18-sqrt(30))/36 (18+sqrt(30))/36 (18+sqrt(30))/36 (18-sqrt(30))/36]/2);
%P=diag([(322-13*sqrt(70))/900 (322+13*sqrt(70))/900 128/225 (322+13*sqrt(70))/900 (322-13*sqrt(70))/900]/2);
P=kron(P,P);
Ohm=kron(P,W);

gridx0=(fullfact(21*ones(1,dx))-1)/20; ngridx0=size(gridx0,1);
gridx2=(fullfact([21 21])-1)/20; ngridx2=size(gridx2,1); 
gridx3=(fullfact(21*ones(1,dx))-1)/20; 
gridz=QuadPoints; ngridz=size(gridz,1);
%Shubert Function + Schwefel Function
targetset=[0.0133379278857657          1.99432496185344;
           13.7858976716857          16.4466194466848];

for i=1:nodesigns
    Designs{i}(:,(dx+1):(dx+dz))=TruncNormInv(Designs{i}(:,(dx+1):(dx+dz)),noiselevel);
    Designs2{i}(:,(dx+1):(dx+dz))=TruncNormInv(Designs2{i}(:,(dx+1):(dx+dz)),noiselevel);  
end

Nrun=size(targetset,2)*nodesigns; 
targetset=kron(targetset,ones(1,nodesigns));

Performance11=cell(1,Nrun); Performance12=cell(1,Nrun); Performance13=cell(1,Nrun); Performance14=cell(1,Nrun); Performance2=cell(1,Nrun);
Distance11=zeros(nadd,Nrun); Distance12=zeros(nadd,Nrun); Distance13=zeros(nadd,Nrun); Distance14=zeros(nadd,Nrun); Distance2=zeros(nadd,Nrun);
for i=1:Nrun
    Performance11{i}=zeros(nadd,4);
    Performance12{i}=zeros(nadd,4);
    Performance13{i}=zeros(nadd,4);      
    Performance14{i}=zeros(nadd,4);          
    Performance2{i}=zeros(1,4);    
end

for run=1:Nrun
target=targetset(:,run); target2=repmat(target,m,1);
trueEQL=zeros(ngridx2,1);
for i=1:ngridx2
    YQ=evalfunc([repmat(gridx2(i,:),m,1) QuadPoints]);
    Dev=YQ-repmat(target,1,m);
    trueEQL(i)=sum(sum(Dev.*(W*Dev*P)));
end
DesignNo=mod(run,nodesigns); 
if(DesignNo==0); DesignNo=nodesigns; end
if(DesignNo==1)
x1=unique(gridx2(:,1)); x2=unique(gridx2(:,2)); 
figure(100+run+(NL-1)*Nrun),mesh(x1,x2,reshape(trueEQL,21,21)');
end
[minEQL,i0]=min(trueEQL); start=gridx2(i0,:);    
[optimalx,optimalEQL,exitflag] = patternsearch(@TrueEQL,start,[],[],[],[],zeros(1,dx),ones(1,dx),[],options);

n=n0;     
D11=Designs{DesignNo};
D12=D11;
D13=D11;
D14=D11;

Y11=evalfunc(D11);
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

while(stop==0)
count=count+1;    

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
parfor i=1:ngridx
    storecriterion11(i)=ComputeEI11(gridx(i,:));
end
[~,i0]=sort(storecriterion11); 
optx=zeros(nostart,dx); fval=zeros(nostart,1);
parfor j=1:nostart
    [optx(j,:),fval(j),exitflag] = patternsearch(ComputeEI11,gridx(i0(j),:),[],[],[],[],zeros(1,dx),ones(1,dx),[],options);
end
[mEI,JJ]=min(fval); optx11=optx(JJ,:);
EI11=[EI11; -mEI]
if((count==No)&&(run==1))
    storecriterion11=zeros(ngridx0,1); 
parfor i=1:ngridx0
    storecriterion11(i)=ComputeEI11(gridx0(i,:));
end    
figure(501),contour([0:0.05:1]',[0:0.05:1]',reshape(storecriterion11(1:ngridx0),21,21)','ShowText','on')
hold on
plot(optx11(1),optx11(2),'^')
input('press enter')
end
[EVar,i0]=CompPV(D11,thetaopt11,alphahat11,sigmahat11,GammaDDI11,Beta11,optx11,target,QuadPoints,m,P,W);
optz11=QuadPoints(i0,:);
D11=[D11; optx11 optz11]
Y11=[Y11 evalfunc(D11(end,:))];

    tau12=minQ12; 
    gridx=[gridx0; D12(:,1:dx); [estoptx; estoptx12(1:(count-1),:)]]; 
    ComputeEI12=@(x)CompEI12(x,D12,thetaopt12,alphahat12,sigmahat12,GammaDDI12,Beta12,target,QuadPoints,P,W,Ohm,m,q,tau12,gridz,ngridz);
    if(count==1)
        Baseline=-ComputeEI12(estoptx);
        Decision=CompareEQL(D12,gridx,QuadPoints,P,W,target,thetaopt12,alphahat12,sigmahat12,GammaDDI12,Beta12,tau12,Baseline);        
    else
        Baseline=-ComputeEI12(estoptx12(count-1,:));
        Decision=CompareEQL(D12,gridx,QuadPoints,P,W,target,thetaopt12,alphahat12,sigmahat12,GammaDDI12,Beta12,tau12,Baseline);
    end    
    gridx=gridx(Decision,:); ngridx=size(gridx,1);    
    storecriterion12=zeros(ngridx,1); 
parfor i=1:ngridx
    storecriterion12(i)=ComputeEI12(gridx(i,:));
end
[~,i0]=sort(storecriterion12);
optx=zeros(nostart,dx); fval=zeros(nostart,1);
for j=1:nostart
    [optx(j,:),fval(j),exitflag] = patternsearch(ComputeEI12,gridx(i0(j),:),[],[],[],[],zeros(1,dx),ones(1,dx),[],options);
end
[mEI,JJ]=min(fval); [mEI,II]=ComputeEI12(optx(JJ,:));
optxz12=[optx(JJ,:) gridz(II,:)];
EI12=[EI12; -mEI]
if((count==No)&&(run==1))
    storecriterion12=zeros(ngridx0,1); 
parfor i=1:ngridx0
    storecriterion12(i)=ComputeEI12(gridx0(i,:));
end        
figure(502),contour([0:0.05:1],[0:0.05:1],reshape(storecriterion12(1:ngridx0),21,21)','ShowText','on')
hold on
plot(optxz12(1),optxz12(2),'^')
input('press enter')
end
D12=[D12; optxz12]
Y12=[Y12 evalfunc(D12(end,:))];

    gridx=[gridx0; D13(:,1:dx); [estoptx; estoptx13(1:(count-1),:)]]; 
    ComputeEI13=@(x)CompEI13(x,D13,thetaopt13,alphahat13,sigmahat13,GammaDDI13,Beta13,target2,QuadPoints,P,W,m); 
    if(count==1)
        Baseline=ComputeEI13(estoptx);
        Decision=CompareEQLV2(D13,gridx,QuadPoints,P,W,target,thetaopt13,alphahat13,sigmahat13,GammaDDI13,Beta13,Baseline);        
    else
        Baseline=ComputeEI13(estoptx13(count-1,:));
        Decision=CompareEQLV2(D13,gridx,QuadPoints,P,W,target,thetaopt13,alphahat13,sigmahat13,GammaDDI13,Beta13,Baseline);
    end    
    gridx=gridx(Decision,:); ngridx=size(gridx,1);    
    storecriterion13=zeros(ngridx,1); 
parfor i=1:ngridx
    storecriterion13(i)=ComputeEI13(gridx(i,:));
end
[~,i0]=sort(storecriterion13); 
optx=zeros(nostart,dx); fval=zeros(nostart,1);
parfor j=1:nostart
    [optx(j,:),fval(j),exitflag] = patternsearch(ComputeEI13,gridx(i0(j),:),[],[],[],[],zeros(1,dx),ones(1,dx),[],options);
end
[mEI,JJ]=min(fval); optx13=optx(JJ,:);
EI13=[EI13; mEI]
if((count==No)&&(run==1))
    storecriterion13=zeros(ngridx0,1); 
parfor i=1:ngridx0
    storecriterion13(i)=ComputeEI13(gridx0(i,:));
end        
figure(503),contour([0:0.05:1]',[0:0.05:1]',reshape(storecriterion13(1:ngridx0),21,21)','ShowText','on')
hold on
plot(optx13(1),optx13(2),'^')
input('press enter')
end
[EVar,i0]=CompPV13(D13,thetaopt13,alphahat13,sigmahat13,GammaDDI13,Beta13,optx13,target,QuadPoints,m,P,W);
optz13=QuadPoints(i0,:);
D13=[D13; optx13 optz13]
Y13=[Y13 evalfunc(D13(end,:))];

    tau14=estobjval14;
    gridx=[gridx0; D14(:,1:dx); [estoptx; estoptx14(1:(count-1),:)]]; 
    ComputeEI14=@(x)CompEI(x,D14,thetaopt14,alphahat14,sigmahat14,GammaDDI14,Beta14,target2,QuadPoints,P,W,m,tau14);   
    if(count==1)
        Baseline=-ComputeEI14(estoptx);
        Decision=CompareEQL(D14,gridx,QuadPoints,P,W,target,thetaopt14,alphahat14,sigmahat14,GammaDDI14,Beta14,tau14,Baseline);        
    else
        Baseline=-ComputeEI14(estoptx14(count-1,:));
        Decision=CompareEQL(D14,gridx,QuadPoints,P,W,target,thetaopt14,alphahat14,sigmahat14,GammaDDI14,Beta14,tau14,Baseline);
    end
    gridx=gridx(Decision,:); ngridx=size(gridx,1);
    storecriterion14=zeros(ngridx,1); 
parfor i=1:ngridx
    storecriterion14(i)=ComputeEI14(gridx(i,:));
end
[~,i0]=sort(storecriterion14); 
optx=zeros(nostart,dx); fval=zeros(nostart,1);
parfor j=1:nostart
    [optx(j,:),fval(j),exitflag] = patternsearch(ComputeEI14,gridx(i0(j),:),[],[],[],[],zeros(1,dx),ones(1,dx),[],options);
end
[mEI,JJ]=min(fval); optx14=optx(JJ,:);
EI14=[EI14; -mEI]
if((count==No)&&(run==1))
    storecriterion14=zeros(ngridx0,1); 
parfor i=1:ngridx0
    storecriterion14(i)=ComputeEI14(gridx0(i,:));
end    
figure(504),contour([0:0.05:1]',[0:0.05:1]',reshape(storecriterion14(1:ngridx0),21,21)','ShowText','on')
hold on
plot(optx14(1),optx14(2),'^')
input('press enter')
end
[EVar,i0]=CompPV(D14,thetaopt14,alphahat14,sigmahat14,GammaDDI14,Beta14,optx14,target,QuadPoints,m,P,W);
optz14=QuadPoints(i0,:);
D14=[D14; optx14 optz14]
Y14=[Y14 evalfunc(D14(end,:))];

n=n+1;

[thetaopt11,alphahat11,sigmahat11,GammaDDI11,Beta11]=GPFitXZV3(D11,Y11,thetaopt11,alphahat11,sigmahat11);
[thetaopt12,alphahat12,sigmahat12,GammaDDI12,Beta12]=GPFitXZV3(D12,Y12,thetaopt12,alphahat12,sigmahat12);
[thetaopt13,alphahat13,sigmahat13,GammaDDI13,Beta13]=GPFitXZV3(D13,Y13,thetaopt13,alphahat13,sigmahat13);
[thetaopt14,alphahat14,sigmahat14,GammaDDI14,Beta14]=GPFitXZV3(D14,Y14,thetaopt14,alphahat14,sigmahat14);

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
if((count==No)&&(run==1))
Objval=zeros(121,1);
gx=(fullfact([11 11])-1)/10;
for i=1:121
    Objval(i)=EQL11(gx(i,:));
end
figure(1001),contour([0:0.1:1]',[0:0.1:1]',reshape(Objval,11,11)','ShowText','on')
hold on
plot(estoptx11(count,1),estoptx11(count,2),'^')
[minObjval,I]=min(Objval);
[gx(I,:) minObjval; estoptx11(count,:) estobjval11]
input('press enter')
end
trueobjval11=TrueEQL(estoptx11(count,:));
estoptx11V2=D11(index,1:dx); estobjval11V2=minQ11; trueobjval11V2=TrueEQL(estoptx11V2);

EQL12=@(x)compEQL(D12,x,QuadPoints,P,W,target,thetaopt12,alphahat12,sigmahat12,GammaDDI12,Beta12);
Q12=EQL12(D12(:,1:dx));
[minQ12,index]=min(Q12);
candidates=[gridx3; estoptx12(1:(count-1),:)];
Q0=EQL12(candidates);
[minQ0,index0]=min(Q0);
if(minQ0<minQ12)
    start=candidates(index0,:);
else
    start=D12(index,1:dx);
end
[estoptx12(count,:),estobjval12,exitflag]=patternsearch(EQL12,start,[],[],[],[],zeros(1,dx),ones(1,dx),[],options); 
if((count==No)&&(run==1))
Objval=zeros(121,1);
gx=(fullfact([11 11])-1)/10;
for i=1:121
    Objval(i)=EQL12(gx(i,:));
end
figure(1002),contour([0:0.1:1]',[0:0.1:1]',reshape(Objval,11,11)','ShowText','on')
hold on
plot(estoptx12(count,1),estoptx12(count,2),'^')
[minObjval,I]=min(Objval);
[gx(I,:) minObjval; estoptx12(count,:) estobjval12]
input('press enter')
end
trueobjval12=TrueEQL(estoptx12(count,:));
estoptx12V2=D12(index,1:dx); estobjval12V2=minQ12; trueobjval12V2=TrueEQL(estoptx12V2);

EQL13=@(x)compEQL(D13,x,QuadPoints,P,W,target,thetaopt13,alphahat13,sigmahat13,GammaDDI13,Beta13);
Q13=EQL13(D13(:,1:dx));
[minQ13,index]=min(Q13);
candidates=[gridx3; estoptx13(1:(count-1),:)];
Q0=EQL13(candidates);
[minQ0,index0]=min(Q0);
if(minQ0<minQ13)
    start=candidates(index0,:);
else
    start=D13(index,1:dx);
end
[estoptx13(count,:),estobjval13,exitflag]=patternsearch(EQL13,start,[],[],[],[],zeros(1,dx),ones(1,dx),[],options); 
if((count==No)&&(run==1))
Objval=zeros(121,1);
gx=(fullfact([11 11])-1)/10;
for i=1:121
    Objval(i)=EQL13(gx(i,:));
end
figure(1003),contour([0:0.1:1]',[0:0.1:1]',reshape(Objval,11,11)','ShowText','on')
hold on
plot(estoptx13(count,1),estoptx13(count,2),'^')
[minObjval,I]=min(Objval);
[gx(I,:) minObjval; estoptx13(count,:) estobjval13]
input('press enter')
end
trueobjval13=TrueEQL(estoptx13(count,:));
estoptx13V2=D13(index,1:dx); estobjval13V2=minQ13; trueobjval13V2=TrueEQL(estoptx13V2);

EQL14=@(x)compEQL(D14,x,QuadPoints,P,W,target,thetaopt14,alphahat14,sigmahat14,GammaDDI14,Beta14);
Q14=EQL14(D14(:,1:dx));
[minQ14,index]=min(Q14);
candidates=[gridx3; estoptx14(1:(count-1),:)];
Q0=EQL14(candidates);
[minQ0,index0]=min(Q0);
if(minQ0<minQ14)
    start=candidates(index0,:);
else
    start=D14(index,1:dx);
end
[estoptx14(count,:),estobjval14,exitflag]=patternsearch(EQL14,start,[],[],[],[],zeros(1,dx),ones(1,dx),[],options); 
if((count==No)&&(run==1))
Objval=zeros(121,1);
gx=(fullfact([11 11])-1)/10;
for i=1:121
    Objval(i)=EQL14(gx(i,:));
end
figure(1004),contour([0:0.1:1]',[0:0.1:1]',reshape(Objval,11,11)','ShowText','on')
hold on
plot(estoptx14(count,1),estoptx14(count,2),'^')
[minObjval,I]=min(Objval);
[gx(I,:) minObjval; estoptx14(count,:) estobjval14]
input('press enter')
end
trueobjval14=TrueEQL(estoptx14(count,:));
estoptx14V2=D14(index,1:dx); estobjval14V2=minQ14; trueobjval14V2=TrueEQL(estoptx14V2);

[(estobjval11-trueobjval11) (trueobjval11-optimalEQL) (estobjval11V2-trueobjval11V2) (trueobjval11V2-optimalEQL)]
[(estobjval12-trueobjval12) (trueobjval12-optimalEQL) (estobjval12V2-trueobjval12V2) (trueobjval12V2-optimalEQL)]
[(estobjval13-trueobjval13) (trueobjval13-optimalEQL) (estobjval13V2-trueobjval13V2) (trueobjval13V2-optimalEQL)]
[(estobjval14-trueobjval14) (trueobjval14-optimalEQL) (estobjval14V2-trueobjval14V2) (trueobjval14V2-optimalEQL)]

Performance11{run}(count,:)=[(estobjval11-trueobjval11) (trueobjval11-optimalEQL) (estobjval11V2-trueobjval11V2) (trueobjval11V2-optimalEQL)];
Performance12{run}(count,:)=[(estobjval12-trueobjval12) (trueobjval12-optimalEQL) (estobjval12V2-trueobjval12V2) (trueobjval12V2-optimalEQL)];
Performance13{run}(count,:)=[(estobjval13-trueobjval13) (trueobjval13-optimalEQL) (estobjval13V2-trueobjval13V2) (trueobjval13V2-optimalEQL)];
Performance14{run}(count,:)=[(estobjval14-trueobjval14) (trueobjval14-optimalEQL) (estobjval14V2-trueobjval14V2) (trueobjval14V2-optimalEQL)];

if(count==nadd)
    stop=1;
end
end

Distance11(:,run)=pdist2(D11((n0+1):end,1:dx),optimalx);
Distance12(:,run)=pdist2(D12((n0+1):end,1:dx),optimalx);
Distance13(:,run)=pdist2(D13((n0+1):end,1:dx),optimalx);
Distance14(:,run)=pdist2(D14((n0+1):end,1:dx),optimalx);

D2=Designs2{DesignNo};
Y2=evalfunc(D2);
[thetaopt2,alphahat2,sigmahat2,GammaDDI2,Beta2]=GPFitXZV3(D2,Y2);
EQL2=@(x)compEQL(D2,x,QuadPoints,P,W,target,thetaopt2,alphahat2,sigmahat2,GammaDDI2,Beta2);
Q2=EQL2(D2(:,1:dx));
[minQ2,index]=min(Q2);
Q0=EQL2(gridx3);
[minQ0,index0]=min(Q0);
if(minQ0<minQ2)
    start=gridx3(index0,:);
else
    start=D2(index,1:dx);
end
[estoptx2,estobjval2,exitflag]=patternsearch(EQL2,start,[],[],[],[],zeros(1,dx),ones(1,dx),[],options);
trueobjval2=TrueEQL(estoptx2);
estoptx2V2=D2(index,1:dx); estobjval2V2=minQ2; trueobjval2V2=TrueEQL(estoptx2V2);
Performance2{run}=[(estobjval2-trueobjval2) (trueobjval2-optimalEQL) (estobjval2V2-trueobjval2V2) (trueobjval2V2-optimalEQL)];
Distance2(:,run)=pdist2(D2((n0+1):end,1:dx),optimalx);

save(strcat('ComparisonSimulationTrial',num2str(run),'noise',num2str(noiselevel),'.mat'))
end

Per11=zeros(nadd,4,Nrun); Per12=zeros(nadd,4,Nrun); Per13=zeros(nadd,4,Nrun); Per14=zeros(nadd,4,Nrun);
for i=1:Nrun
    Per11(:,:,i)=Performance11{i}; Per12(:,:,i)=Performance12{i}; Per13(:,:,i)=Performance13{i}; Per14(:,:,i)=Performance14{i};    
end

figure(1+3*(NL-1)),plot(1:nadd,mean(Per11(:,2,1:nodesigns),3))
hold on,
plot(1:nadd,mean(Per12(:,2,1:nodesigns),3))
plot(1:nadd,mean(Per13(:,2,1:nodesigns),3))
plot(1:nadd,mean(Per14(:,2,1:nodesigns),3))
legend('1','2','3','4')

figure(2+3*(NL-1)),plot(1:nadd,mean(Per11(:,2,(nodesigns+1):(2*nodesigns)),3))
hold on,
plot(1:nadd,mean(Per12(:,2,(nodesigns+1):(2*nodesigns)),3))
plot(1:nadd,mean(Per13(:,2,(nodesigns+1):(2*nodesigns)),3))
plot(1:nadd,mean(Per14(:,2,(nodesigns+1):(2*nodesigns)),3))
legend('1','2','3','4')

end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function EQL=TrueEQL(x)
global target QuadPoints P W m
YQ=evalfunc([repmat(x,m,1) QuadPoints]);
Dev=YQ-repmat(target,1,m);
EQL=sum(sum(Dev.*(W*Dev*P)));
end

function [mEI]=CompEI(x,D,thetaopt,alphahat,sigmahat,GammaDDI,Beta,target2,QuadPoints,P,W,m,tau)
XZ=[kron(x,ones(m,1)) QuadPoints];
[M,V2]=GPPredictXZ(D,thetaopt,alphahat,sigmahat,GammaDDI,Beta,XZ);
M=M(:); 
EI=EIQ(M,V2,sigmahat,P,W,target2,tau);
mEI=-EI;
end

function [EVar,II]=CompPV(D,thetaopt,alphahat,sigmahat,GammaDDI,Beta,optx,target,QuadPoints,m,P,W)

[M,V2]=GPPredictXZ(D,thetaopt,alphahat,sigmahat,GammaDDI,Beta,[repmat(optx,m,1) QuadPoints]);
EVar=zeros(m,1);
parfor i=1:m
    Vi=V2; Ri=Vi(:,i); 
    CVi=Vi-Ri*(Ri')/(Vi(i,i)+10^-6);
    CVi(i,:)=zeros(1,m); CVi(:,i)=zeros(m,1);
    Dev=M-repmat(target,1,m);    
    A=P*CVi; B=W*sigmahat; C=B*W*Dev*(A*P);
    EVar(i)=2*sum(sum(A.*A,1))*sum(sum(B.*B,1))+4*Dev(:)'*C(:)+4*(Ri'/Vi(i,i)*A*P*Ri/Vi(i,i))*sum(sum((B*W).*sigmahat,1))*Vi(i,i);
end
[minEVar,II]=min(EVar);
end

function [mEI]=CompEI13(x,D13,thetaopt13,alphahat13,sigmahat13,GammaDDI13,Beta13,target2,QuadPoints,P,W,m)
XZ=[kron(x,ones(m,1)) QuadPoints];
[M,V2]=GPPredictXZ(D13,thetaopt13,alphahat13,sigmahat13,GammaDDI13,Beta13,XZ);
M=M(:); 
mEI=ComputeQuantile(M,V2,sigmahat13,P,W,target2);
end

function [EVar,II]=CompPV13(D13,thetaopt13,alphahat13,sigmahat13,GammaDDI13,Beta13,optx13,target,QuadPoints,m,P,W)

[M,V2]=GPPredictXZ(D13,thetaopt13,alphahat13,sigmahat13,GammaDDI13,Beta13,[repmat(optx13,m,1) QuadPoints]);
EVar=zeros(m,1);
parfor i=1:m
    Vi=V2; Ri=Vi(:,i); 
    CVi=Vi-Ri*(Ri')/(Vi(i,i)+10^-6);
    CVi(i,:)=zeros(1,m); CVi(:,i)=zeros(m,1);
    Dev=M-repmat(target,1,m);    
    A=P*CVi; B=W*sigmahat13; C=B*W*Dev*(A*P);
    EVar(i)=2*sum(sum(A.*A,1))*sum(sum(B.*B,1))+4*Dev(:)'*C(:)+4*(Ri'/Vi(i,i)*A*P*Ri/Vi(i,i))*sum(sum((B*W).*sigmahat13,1))*Vi(i,i);
end
[minEVar,II]=min(EVar);
end

function [minmEI,II]=CompEI12(x,D12,thetaopt12,alphahat12,sigmahat12,GammaDDI12,Beta12,target,QuadPoints,P,W,Ohm,m,q,tau12,gridz,ngridz)
mEI=zeros(ngridz,1);
xz=[repmat(x,ngridz,1) gridz];
[M3,V3]=GPPredictXZ(D12,thetaopt12,alphahat12,sigmahat12,GammaDDI12,Beta12,[repmat(x,m,1) QuadPoints; xz]);
parfor jj=1:ngridz
M2=M3; V2=V3;
M2=M2(:,[1:m,m+jj]); V2=V2([1:m,m+jj],[1:m,m+jj]);
V2(m+1,m+1)=V2(m+1,m+1)+10^-6;
astar=zeros(q,m); 
b=zeros(m,1); 
newPsi=zeros(m,1); 
for i=1:m
    rho=V2(i,m+1)/V2(m+1,m+1);
    astar(:,i)=M2(:,i)-target-rho*M2(:,m+1);
    b(i)=+rho;
    newPsi(i)=max(V2(i,i)-V2(i,m+1)^2/V2(m+1,m+1),0); 
end

av=astar(:); bTPb=b'*P*b; c=astar*P*b; e=c/bTPb;
l=-tau12+av'*Ohm*av-c'*W*c/bTPb+sum(diag(P).*newPsi)*trace(W*sigmahat12);
if((sum(b~=0)>0)&&(l<0))
    M=M2(:,m+1)+e;
    C=V2(m+1,m+1)*sigmahat12;
    E=bTPb*W;
    mEI(jj)=-EMIQ(M,C,E,l);
elseif(l>=0)
    mEI(jj)=0;
else
    l=-tau12+av'*Ohm*av+sum(diag(P).*newPsi)*trace(W*sigmahat12);
    mEI(jj)=-max(-l,0);
end
end
[minmEI,II]=min(mEI);

end

function Q=compEQL(D,x,QuadPoints,P,W,target,thetaopt,alphahat,sigmahat,GammaDDI,Beta)
m=size(QuadPoints,1);
N=size(x,1);
Q=zeros(N,1);
parfor i=1:N
    XZ=[kron(x(i,:),ones(m,1)) QuadPoints];
    [M,V2]=GPPredictXZ(D,thetaopt,alphahat,sigmahat,GammaDDI,Beta,XZ);
    Dev=M-repmat(target,1,m);
    Q(i)=sum(sum(Dev.*(W*Dev*P)))+sum(diag(P).*diag(V2))*sum(sum(W.*sigmahat,1));
end
end

function Decision=CompareEQL(D,x,QuadPoints,P,W,target,thetaopt,alphahat,sigmahat,GammaDDI,Beta,tau,Baseline)
m=size(QuadPoints,1);
N=size(x,1);
Decision=zeros(N,1);
parfor i=1:N
    XZ=[kron(x(i,:),ones(m,1)) QuadPoints];
    [M,V2]=GPPredictXZ(D,thetaopt,alphahat,sigmahat,GammaDDI,Beta,XZ);
    Dev=M-repmat(target,1,m);    
    Mean=sum(sum(Dev.*(W*Dev*P)))+sum(diag(P).*diag(V2))*sum(sum(W.*sigmahat,1));
    A=P*V2; B=W*sigmahat; C=B*W*Dev*(A*P);
    Std=sqrt(max(2*sum(sum(A.*A,1))*sum(sum(B.*B,1))+4*Dev(:)'*C(:),0));
    k=max([(Mean-tau)/Std,0]);
    if(k>0)
        UpperBound=tau/(1+k^2);
    else
        UpperBound=Baseline+1;
    end
    if(UpperBound<Baseline)
        Decision(i)=0;
    else
        Decision(i)=1;
    end
end
    Decision=logical(Decision);
end

function Decision=CompareEQLV2(D,x,QuadPoints,P,W,target,thetaopt,alphahat,sigmahat,GammaDDI,Beta,Baseline)
m=size(QuadPoints,1);
N=size(x,1);
Decision=zeros(N,1);
parfor i=1:N
    XZ=[kron(x(i,:),ones(m,1)) QuadPoints];
    [M,V2]=GPPredictXZ(D,thetaopt,alphahat,sigmahat,GammaDDI,Beta,XZ);
    Dev=M-repmat(target,1,m);    
    Mean=sum(sum(Dev.*(W*Dev*P)))+sum(diag(P).*diag(V2))*sum(sum(W.*sigmahat,1));
    A=P*V2; B=W*sigmahat; C=B*W*Dev*(A*P);
    Std=sqrt(max(2*sum(sum(A.*A,1))*sum(sum(B.*B,1))+4*Dev(:)'*C(:),0));
    Thres=Mean-sqrt(1/0.025-1)*Std;
    if(Thres>Baseline)
        Decision(i)=0;
    else
        Decision(i)=1;
    end
end
    Decision=logical(Decision);
end

function Y=evalfunc(D)
%Shubert Function + Schwefel Function (*)
u1=1*(D(:,1)+D(:,3))+1; v1=1*(D(:,2)+D(:,4))+1;
n=size(u1,1); 
u2=50*(D(:,1)+D(:,3))+50; v2=50*(D(:,2)+D(:,4))+50;
Y=[(sum(repmat(1:5,n,1).*cos(repmat((1:5)+1,n,1).*repmat(u1,1,5)+repmat((1:5),n,1)),2).*sum(repmat(1:5,n,1).*cos(repmat((1:5)+1,n,1).*repmat(v1,1,5)+repmat((1:5),n,1)),2))';
    (418.9829*2-u2.*sin(sqrt(abs(u2)))-v2.*sin(sqrt(abs(v2))))']./[3.72858741286651;57.5662387791625];
end

function ComparisonSimulationBayesOptimizationd2ParallelIR
%checked$27042019
clear all
format long g
global dx target target2 QuadPoints P W Ohm m q 
%checked$
No=Inf; 
noiselevels=[0.15]; lnoiselevels=length(noiselevels); nostart=1; options=optimoptions(@patternsearch,'MaxIter',10000,'TolFun',10^-9);
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
%Shubert Function + Schwefel Function
targetset=[0.0133379278857657          1.99432496185344;
           13.7858976716857          16.4466194466848];

for i=1:nodesigns
    Designs{i}(:,(dx+1):(dx+dz))=TruncNormInv(Designs{i}(:,(dx+1):(dx+dz)),noiselevel);
end

Nrun=size(targetset,2)*nodesigns; 
targetset=kron(targetset,ones(1,nodesigns));

Performance15=cell(1,Nrun); Distance15=zeros(nadd,Nrun); 
for i=1:Nrun
    Performance15{i}=zeros(nadd,4);
end

for run=1:Nrun
target=targetset(:,run); target2=repmat(target,m,1);
trueEQL=zeros(ngridx2,1);
for i=1:ngridx2
    YQ=evalfunc([repmat(gridx2(i,:),m,1) QuadPoints]);
    trueEQL(i)=sum(diag(P).*YQ(:));
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
D15=Designs{DesignNo};

Y15=evalfunc(D15);

stop=0; count=0; EI15=[]; estoptx15=zeros(nadd,dx); 

[thetaopt15,alphahat15,sigmahat15,GammaDDI15,Beta15]=GPFitXZV3(D15,Y15);

Q15=compEQL(D15,D15(:,1:dx),QuadPoints,P,thetaopt15,alphahat15,sigmahat15,GammaDDI15,Beta15);
[minQ15,index]=min(Q15);
Q0=compEQL(D15,gridx3,QuadPoints,P,thetaopt15,alphahat15,sigmahat15,GammaDDI15,Beta15);
[minQ0,index0]=min(Q0);
if(minQ0<minQ15)
    start=gridx3(index0,:);
else
    start=D15(index,1:dx);
end
EQL15=@(x)compEQL(D15,x,QuadPoints,P,thetaopt15,alphahat15,sigmahat15,GammaDDI15,Beta15);        
[estoptx,estobjval15,exitflag]=patternsearch(EQL15,start,[],[],[],[],zeros(1,dx),ones(1,dx),[],options);         

while(stop==0)
count=count+1;    

    tau15=minQ15; 
    gridx=[gridx0; D15(:,1:dx); [estoptx; estoptx15(1:(count-1),:)]]; ngridx=size(gridx,1);
    ComputeEI15=@(x)CompEI(x,D15,thetaopt15,alphahat15,sigmahat15,GammaDDI15,Beta15,QuadPoints,P,m,tau15);   
    storecriterion15=zeros(ngridx,1); 
parfor i=1:ngridx
    storecriterion15(i)=ComputeEI15(gridx(i,:));
end
[~,i0]=sort(storecriterion15); 
optx=zeros(nostart,dx); fval=zeros(nostart,1);
parfor j=1:nostart
    [optx(j,:),fval(j),exitflag] = patternsearch(ComputeEI15,gridx(i0(j),:),[],[],[],[],zeros(1,dx),ones(1,dx),[],options);
end
[mEI,JJ]=min(fval); optx15=optx(JJ,:);
EI15=[EI15; -mEI]
if((count==No)&&(run==1))
    storecriterion15=zeros(ngridx0,1); 
parfor i=1:ngridx0
    storecriterion15(i)=ComputeEI15(gridx0(i,:));
end    
figure(501),contour([0:0.05:1]',[0:0.05:1]',reshape(storecriterion15(1:ngridx0),21,21)','ShowText','on')
hold on
plot(optx15(1),optx15(2),'^')
input('press enter')
end
[EVar,i0]=CompPV(D15,thetaopt15,alphahat15,sigmahat15,GammaDDI15,Beta15,optx15,QuadPoints,m,P);
optz15=QuadPoints(i0,:);
D15=[D15; optx15 optz15]
Y15=[Y15 evalfunc(D15(end,:))];

n=n+1;

[thetaopt15,alphahat15,sigmahat15,GammaDDI15,Beta15]=GPFitXZV3(D15,Y15,thetaopt15,alphahat15,sigmahat15);

EQL15=@(x)compEQL(D15,x,QuadPoints,P,thetaopt15,alphahat15,sigmahat15,GammaDDI15,Beta15);
Q15=EQL15(D15(:,1:dx));
[minQ15,index]=min(Q15);
candidates=[gridx3; estoptx15(1:(count-1),:)];
Q0=EQL15(candidates);
[minQ0,index0]=min(Q0);
if(minQ0<minQ15)
    start=candidates(index0,:);
else
    start=D15(index,1:dx);
end
[estoptx15(count,:),estobjval15,exitflag]=patternsearch(EQL15,start,[],[],[],[],zeros(1,dx),ones(1,dx),[],options); 
if((count==No)&&(run==1))
Objval=zeros(121,1);
gx=(fullfact([11 11])-1)/10;
for i=1:121
    Objval(i)=EQL15(gx(i,:));
end
figure(1001),contour([0:0.1:1]',[0:0.1:1]',reshape(Objval,11,11)','ShowText','on')
hold on
plot(estoptx15(count,1),estoptx15(count,2),'^')
[minObjval,I]=min(Objval);
[gx(I,:) minObjval; estoptx15(count,:) estobjval15]
input('press enter')
end
trueobjval15=TrueEQL(estoptx15(count,:));
estoptx15V2=D15(index,1:dx); estobjval15V2=minQ15; trueobjval15V2=TrueEQL(estoptx15V2);

[(estobjval15-trueobjval15) (trueobjval15-optimalEQL) (estobjval15V2-trueobjval15V2) (trueobjval15V2-optimalEQL)]

Performance15{run}(count,:)=[(estobjval15-trueobjval15) (trueobjval15-optimalEQL) (estobjval15V2-trueobjval15V2) (trueobjval15V2-optimalEQL)];

if(count==nadd)
    stop=1;
end
end

Distance15(:,run)=pdist2(D15((n0+1):end,1:dx),optimalx);

save(strcat('IRSimulationTrial',num2str(run),'noise',num2str(noiselevel),'.mat'))
end

Per15=zeros(nadd,4,Nrun); 
for i=1:Nrun
    Per15(:,:,i)=Performance15{i}; 
end

figure(1+3*(NL-1)),plot(1:nadd,mean(Per15(:,2,1:nodesigns),3))

figure(2+3*(NL-1)),plot(1:nadd,mean(Per15(:,2,(nodesigns+1):(2*nodesigns)),3))

end
end

function EQL=TrueEQL(x)
global QuadPoints P m
YQ=evalfunc([repmat(x,m,1) QuadPoints]);
EQL=sum(diag(P).*YQ(:));
end

function [mEI]=CompEI(x,D,thetaopt,alphahat,sigmahat,GammaDDI,Beta,QuadPoints,P,m,tau)
XZ=[kron(x,ones(m,1)) QuadPoints];
[M,V2,V]=GPPredictXZ(D,thetaopt,alphahat,sigmahat,GammaDDI,Beta,XZ);
P2=diag(P);
Mu=sum(P2.*M(:));
Std=sqrt(max(P2'*V*P2,0));
if(Std>0)
    mEI=-((tau-Mu)*normcdf((tau-Mu)/Std)+Std*normpdf((tau-Mu)/Std));
else
    mEI=0;
end
end

function [EVar,II]=CompPV(D,thetaopt,alphahat,sigmahat,GammaDDI,Beta,optx,QuadPoints,m,P)

[M,V2,V]=GPPredictXZ(D,thetaopt,alphahat,sigmahat,GammaDDI,Beta,[repmat(optx,m,1) QuadPoints]);
EVar=zeros(m,1); P2=diag(P);
parfor i=1:m
    Vi=V; Ri=Vi(:,i); 
    CVi=Vi-Ri*(Ri')/(Vi(i,i)+10^-6);
    CVi(i,:)=zeros(1,m); CVi(:,i)=zeros(m,1);
    EVar(i)=P2'*CVi*P2;
end
[minEVar,II]=min(EVar);
end

function Q=compEQL(D,x,QuadPoints,P,thetaopt,alphahat,sigmahat,GammaDDI,Beta)
m=size(QuadPoints,1);
N=size(x,1);
Q=zeros(N,1);
parfor i=1:N
    XZ=[kron(x(i,:),ones(m,1)) QuadPoints];
    [M,V2]=GPPredictXZ(D,thetaopt,alphahat,sigmahat,GammaDDI,Beta,XZ);
    Q(i)=sum(diag(P).*M(:));
end
end

function Y=evalfunc(D)
global W target
%Shubert Function + Schwefel Function (*)
u1=1*(D(:,1)+D(:,3))+1; v1=1*(D(:,2)+D(:,4))+1;
n=size(u1,1); 
u2=50*(D(:,1)+D(:,3))+50; v2=50*(D(:,2)+D(:,4))+50;
Y2=[(sum(repmat(1:5,n,1).*cos(repmat((1:5)+1,n,1).*repmat(u1,1,5)+repmat((1:5),n,1)),2).*sum(repmat(1:5,n,1).*cos(repmat((1:5)+1,n,1).*repmat(v1,1,5)+repmat((1:5),n,1)),2))';
    (418.9829*2-u2.*sin(sqrt(abs(u2)))-v2.*sin(sqrt(abs(v2))))']./[3.72858741286651;57.5662387791625];
Dev=Y2-repmat(target,1,n);
Y=sum(Dev.*(W*Dev),1);
end

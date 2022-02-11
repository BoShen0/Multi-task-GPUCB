function [thetaopt,alphahat,sigmahat,GammaDDI,Beta,GammaDD]=GPFitXZV3(Din,Yin,thetaopt0,alphahat0,sigmahat0)
%checked$20052019

global Dfit Yfit nfit dfit qfit scaling C
d=size(Din,2); d=d/2;
Dfit=Din(:,1:d)+Din(:,(d+1):(2*d));
[nfit,dfit]=size(Dfit); 
Yfit=Yin; qfit=size(Yfit,1);
scaling=[range(Dfit,1)];
C=cell(1,dfit);
for i=1:dfit
    C{i}=abs(repmat(Dfit(:,i),1,nfit)-repmat(Dfit(:,i)',nfit,1));
end
options=optimoptions(@patternsearch,'Display','off');
%if(nargin<3)
%[paropt,fval,exitflag]=patternsearch(@Obj,[log((1.44269504088896.*(scaling.^2)).*ones(1,dfit))],[],[],[],[],[log(0.072382413650542.*scaling.^2)],[log(499.75.*scaling.^2)],[],options); 
%else
%[paropt,fval,exitflag]=patternsearch(@Obj,[log(start)],[],[],[],[],[log(0.072382413650542.*scaling.^2)],[log(499.75.*scaling.^2)],[],options); 
%end
if(nargin<3)
    if(d<=3) 
        t=6;
    else
        t=3;
    end
paropt=zeros(t^dfit,dfit); fval=zeros(t^dfit,1);
grid=(fullfact(t*ones(1,dfit))-1)/(t-1); grid=grid*(22.0254528130684-0.108302309059412)+0.108302309059412;
for j=1:(t^dfit)
    [paropt(j,:),fval(j),exitflag]=patternsearch(@Obj,[log((grid(j,:).*scaling).*ones(1,dfit))],[],[],[],[],[log(0.108302309059412.*scaling)],[log(22.0254528130684.*scaling)],[],options); 
end    
[minfval,jj]=min(fval);
thetaopt=[exp(paropt(jj,:))];

disp('Optimal values of GP parameters')
disp(thetaopt)

GammaDD=CompCorr1(thetaopt);
[GammaDDI]=invandlogdet(GammaDD);
GammaDDIOnes=GammaDDI*ones(nfit,1);
OnesGammaDDIOnes=(ones(1,nfit)*GammaDDIOnes);
alphahat=Yfit*GammaDDIOnes/OnesGammaDDIOnes;
Res=Yfit-alphahat*ones(1,nfit);
Beta=Res*GammaDDI;
sigmahat=(Beta*(Res'))/(nfit);

elseif(nargin==3)
    
thetaopt=thetaopt0;
GammaDD=CompCorr1(thetaopt);
[GammaDDI]=invandlogdet(GammaDD);
GammaDDIOnes=GammaDDI*ones(nfit,1);
OnesGammaDDIOnes=(ones(1,nfit)*GammaDDIOnes);
alphahat=Yfit*GammaDDIOnes/OnesGammaDDIOnes;
Res=Yfit-alphahat*ones(1,nfit);
Beta=Res*GammaDDI;
sigmahat=(Beta*(Res'))/(nfit);

else
    
thetaopt=thetaopt0;
alphahat=alphahat0;
sigmahat=sigmahat0;    
GammaDD=CompCorr1(thetaopt);
[GammaDDI]=invandlogdet(GammaDD);
Res=Yfit-alphahat*ones(1,nfit);
Beta=Res*GammaDDI;

end

function Objective=Obj(par)
global Yfit nfit qfit scaling

if(sum(isnan(par))>0)
    Objective=Inf;
    return
end
theta=exp(par);
%if(sum(theta>(499.75.*(scaling.^2)))>0)
%    Objective=Inf;
%    return
%end
if(sum(theta>(22.0254528130684.*scaling))>0)
    Objective=Inf;
    return
end

GammaDD=CompCorr1(theta);
[GammaDDI,LDGammaDD]=invandlogdet(GammaDD);
GammaDDIOnes=GammaDDI*ones(nfit,1);
OnesGammaDDIOnes=(ones(1,nfit)*GammaDDIOnes);
alphahat=Yfit*GammaDDIOnes/OnesGammaDDIOnes;
Res=Yfit-alphahat*ones(1,nfit);
Beta=Res*GammaDDI;
sigmahat=(Beta*(Res'))/(nfit);
Objective=qfit*LDGammaDD+nfit*log(det(sigmahat));

function Corr=CompCorr1(theta)
global nfit dfit C
Corr=1;
for i=1:dfit
%    Corr=Corr.*(exp(-0.5*(C{i}.^2)./theta(i)));    
    rho=C{i}/theta(i);
    Corr=Corr.*(exp(-rho).*(1+rho));
end      
Corr=Corr+10^-6*speye(nfit);
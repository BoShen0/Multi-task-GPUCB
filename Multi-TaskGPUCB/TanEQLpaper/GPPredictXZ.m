function [M,V2,V]=GPPredictXZ(Din,thetaopt,alphahat,sigmahat,GammaDDI,Beta,xz)
%checked$20052019

d=size(Din,2); d=d/2;
D=Din(:,1:d)+Din(:,(d+1):(2*d));
n=size(D,1); 
q=size(alphahat,1); 
L=size(xz,1);
x=xz(:,1:d)+xz(:,(d+1):(2*d));

GammaDx=zeros(n,L); V0=diag(ones(1,L));
for i=1:L
    GammaDx(:,i)=CompCorr1(D,repmat(x(i,:),n,1),thetaopt);
    for j=(i+1):L
        V0(i,j)=CompCorr1(x(i,:),x(j,:),thetaopt);
        V0(j,i)=V0(i,j);
    end
end
M=repmat(alphahat,1,L)+Beta*GammaDx;
V2=V0-GammaDx'*GammaDDI*GammaDx;
V2=(V2+V2')/2;
if(nargout==3)
V=kron(V2,sigmahat);
end

function Corr=CompCorr1(xz1,xz2,theta)
N=size(xz1,1);
%Corr=prod(exp(-0.5*((xz1-xz2).^2)./repmat(theta,N,1)),2);
rho=abs(xz1-xz2)./repmat(theta,N,1);
Corr=prod(exp(-rho).*(1+rho),2);
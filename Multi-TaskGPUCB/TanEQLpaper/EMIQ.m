function EI=EMIQ(M,C,Ohm,lin)
%checked$24052019
global l lambda mu EQ2
warning off all

l=lin;
[eigvecOhm,eigvalOhm]=svd(Ohm); srOhm=eigvecOhm*sqrt(eigvalOhm)*eigvecOhm';
M=srOhm*M; U=srOhm*C*srOhm;

[V,E]=svd(U);
lambda=diag(E); lambda=lambda(:)';
mu=(M)'*V;
Mean=-imag(DCF(0));
r=length(mu);
EQ2=sum(mu.^4+6*lambda.*mu.^2+3*lambda.^2)+(mu.^2+lambda)*ones(r,r)*(mu.^2+lambda)'-sum((mu.^2+lambda).^2)+l^2+2*l*sum(mu.^2+lambda);
StdQ=sqrt(EQ2-Mean^2);
if((abs(Mean)/StdQ)>sqrt(1/0.005-1)) %Chebyshev 
EI=max(Mean,0);
else
EI=Mean/2-(1/(2*pi))*real(integral(@integrand,0,Inf,'AbsTol',1e-9,'RelTol',1e-6));
end

function vals=integrand(t)
global EQ2
t=t(:)';
[DCF1,DCF2]=DCF(t);
vals=(-DCF2+DCF1)./t;
vals(t==0)=-2*(EQ2);

function [vals,vals2]=DCF(t)
global l lambda mu

t=t(:); lt=length(t);
r=length(mu);
t2=repmat(t,1,r);
lambda2=repmat(lambda,lt,1);
mu2=repmat(mu,lt,1);
CF=exp(1i*l*t).*prod((1-2*1i*lambda2.*t2).^(-0.5).*exp((1i*mu2.^2.*t2)./(1-2*1i*lambda2.*t2)),2);
vals=(sum((2*lambda2.^2.*t2+lambda2*1i+(mu2.^2)*1i)./((1-2*1i*lambda2.*t2).^2),2)+1i*l).*CF;
vals=reshape(vals,1,lt);
t2=-t2;
vals2=(sum((2*lambda2.^2.*t2+lambda2*1i+(mu2.^2)*1i)./((1-2*1i*lambda2.*t2).^2),2)+1i*l).*conj(CF);
vals2=reshape(vals2,1,lt);
function EI=EIQ(M,Psi,sigmahat,P,W,target,tauin,boAM)
%checked$04052019
global tau lambda mu 
warning off all

tau=tauin;
[Vm,Em]=svd(P);
srP=Vm*sqrt(Em)*Vm';
Um=srP*Psi*srP;
[gammam,lambdam]=svd(Um);
[Vq,Eq]=svd(W);
srW=Vq*sqrt(Eq)*Vq';
Uq=srW*sigmahat*srW;
[gammaq,lambdaq]=svd(Uq);
lambda=kron(diag(lambdam),diag(lambdaq)); lambda=lambda(:)';
mu=(M-target)'*kron(srP*gammam,srW*gammaq);
if boAM>3
EI=max(real(integral(@Integrand,0,Inf,'AbsTol',1e-6,'RelTol',1e-3)),0);
else
EI=max(real(integral(@Integrand,0,Inf,'AbsTol',1e-9,'RelTol',1e-6)),0);   
% EI=max(real(integral(@Integrand,0,Inf,'AbsTol',1e-9,'RelTol',1e-6)),0);
% it was used for all the two-dimensional functions in the numerical study
end

function vals=Integrand(t)
global tau lambda mu

t=t(:); lt=length(t);
r=length(mu);
t2=repmat(t,1,r);
lambda2=repmat(lambda,lt,1);
mu2=repmat(mu,lt,1);
vals=1/(pi)*((1-exp(-1i*t*tau)-1i*tau*t)./t.^2).*prod((1-2*1i*lambda2.*t2).^(-0.5).*exp((1i*mu2.^2.*t2)./(1-2*1i*lambda2.*t2)),2);
vals(find(t==0))=tau^2/(2*pi);
vals=reshape(vals,1,lt);
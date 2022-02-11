function Quantile=ComputeQuantile(M,Psi,sigmahat,P,W,target)
%checked$04052019
warning off all

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
m=size(P,1); q=size(W,1);
Dev=reshape(M-target,q,m);
x0=sum(sum(Dev.*(W*Dev*P)))+sum(diag(P).*diag(Psi))*sum(sum((W.*sigmahat),1));
vals=zeros(101,1);
r=m*q; Std=sqrt((sum(mu.^4+6*lambda.*mu.^2+3*lambda.^2)+(mu.^2+lambda)*ones(r,r)*(mu.^2+lambda)'-sum((mu.^2+lambda).^2))-x0^2);
l=max(0,x0-sqrt(1/0.025-1)*Std); %Chebyshev 
u=min(x0/0.975,x0+sqrt(1/0.975-1)*Std); %Chebyshev and Markov
grid=[l:((u-l)/100):u]';
DevFunc=@(x)Deviation(x,lambda,mu,l,u);
parfor i=1:101
vals(i)=DevFunc(grid(i));
end
[minvals,I]=min(vals); options=optimoptions(@patternsearch,'Display','off');
[q,fval]=patternsearch(DevFunc,grid(I),[],[],[],[],l,u,[],options);
if(fval>10^-3)
Quantile=grid(I);
else
Quantile=q;
end

function Dev=Deviation(xin,lambda,mu,l,u)
x=xin;
if((x>u)||(x<l))
Dev=1;
else
CFFunc=@(t)CF(t,x,lambda,mu);    
CP=2*real(integral(CFFunc,0,Inf,'AbsTol',1e-9,'RelTol',1e-6));
CP=max(CP,0);
CP=min(CP,1);
Dev=(CP-0.025).^2;
end

function vals=CF(t,x,lambda,mu)
h=x;
t=t(:); lt=length(t);
r=length(mu);
t2=repmat(t,1,r);
lambda2=repmat(lambda,lt,1);
mu2=repmat(mu,lt,1);
vals=1/(pi)*((sin(h*t/2)./t).*exp(-1i*t*h/2)).*prod((1-2*1i*lambda2.*t2).^(-0.5).*exp((1i*mu2.^2.*t2)./(1-2*1i*lambda2.*t2)),2);
vals(find(t==0))=h/(2*pi);
vals=reshape(vals,1,lt);
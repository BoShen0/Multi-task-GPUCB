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

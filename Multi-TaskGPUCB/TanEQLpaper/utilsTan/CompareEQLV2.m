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

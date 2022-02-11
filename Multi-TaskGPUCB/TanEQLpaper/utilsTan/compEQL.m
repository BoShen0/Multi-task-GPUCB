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

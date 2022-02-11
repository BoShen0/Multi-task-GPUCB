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


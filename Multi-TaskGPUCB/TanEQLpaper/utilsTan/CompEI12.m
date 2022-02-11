function [minmEI,II]=CompEI12(x,D12,thetaopt12,alphahat12,sigmahat12,GammaDDI12,Beta12,target,QuadPoints,P,W,Ohm,m,q,tau12,gridz,ngridz)
mEI=zeros(ngridz,1);
xz=[repmat(x,ngridz,1) gridz];
[M3,V3]=GPPredictXZ(D12,thetaopt12,alphahat12,sigmahat12,GammaDDI12,Beta12,[repmat(x,m,1) QuadPoints; xz]);
parfor jj=1:ngridz
M2=M3; V2=V3;
M2=M2(:,[1:m,m+jj]); V2=V2([1:m,m+jj],[1:m,m+jj]);
V2(m+1,m+1)=V2(m+1,m+1)+10^-6;
astar=zeros(q,m); 
b=zeros(m,1); 
newPsi=zeros(m,1); 
for i=1:m
    rho=V2(i,m+1)/V2(m+1,m+1);
    astar(:,i)=M2(:,i)-target-rho*M2(:,m+1);
    b(i)=+rho;
    newPsi(i)=max(V2(i,i)-V2(i,m+1)^2/V2(m+1,m+1),0); 
end

av=astar(:); bTPb=b'*P*b; c=astar*P*b; e=c/bTPb;
l=-tau12+av'*Ohm*av-c'*W*c/bTPb+sum(diag(P).*newPsi)*trace(W*sigmahat12);
if((sum(b~=0)>0)&&(l<0))
    M=M2(:,m+1)+e;
    C=V2(m+1,m+1)*sigmahat12;
    E=bTPb*W;
    mEI(jj)=-EMIQ(M,C,E,l);
elseif(l>=0)
    mEI(jj)=0;
else
    l=-tau12+av'*Ohm*av+sum(diag(P).*newPsi)*trace(W*sigmahat12);
    mEI(jj)=-max(-l,0);
end
end
[minmEI,II]=min(mEI);

end

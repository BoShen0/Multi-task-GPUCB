function [mEI]=CompEI(x,D,thetaopt,alphahat,sigmahat,GammaDDI,Beta,target2,QuadPoints,P,W,m,tau)
% % global m2
% m2 = m;
XZ=[kron(x,ones(m,1)) QuadPoints];
[M,V2]=GPPredictXZ(D,thetaopt,alphahat,sigmahat,GammaDDI,Beta,XZ);
M=M(:); 
EI=EIQ(M,V2,sigmahat,P,W,target2,tau,size(QuadPoints,2));
mEI=-EI;
end


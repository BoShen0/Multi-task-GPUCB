function [mEI]=CompEI13(x,D13,thetaopt13,alphahat13,sigmahat13,GammaDDI13,Beta13,target2,QuadPoints,P,W,m)
XZ=[kron(x,ones(m,1)) QuadPoints];
[M,V2]=GPPredictXZ(D13,thetaopt13,alphahat13,sigmahat13,GammaDDI13,Beta13,XZ);
M=M(:); 
mEI=ComputeQuantile(M,V2,sigmahat13,P,W,target2);
end


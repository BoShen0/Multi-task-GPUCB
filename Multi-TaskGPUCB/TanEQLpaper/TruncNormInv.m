function XT=TruncNormInv(X,Std)

XT=norminv(X*(normcdf(3,0,1)-normcdf(-3,0,1))+normcdf(-3,0,1),0,Std);
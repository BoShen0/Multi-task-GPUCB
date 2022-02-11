function EQL=TrueEQL(x)
global target QuadPoints P W m
YQ=evalfunc([repmat(x,m,1) QuadPoints]);
Dev=YQ-repmat(target,1,m);
EQL=sum(sum(Dev.*(W*Dev*P)));
end


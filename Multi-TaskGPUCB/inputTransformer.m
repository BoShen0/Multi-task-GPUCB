function xnew = inputTransformer(x, boundsSource, boundsTarget)
% this function is used to transform
Ratio = (x' - boundsSource(:,1))./(boundsSource(:,2) - boundsSource(:,1));
 xprime = boundsTarget(:,1) + (boundsTarget(:,2) - boundsTarget(:,1)).* Ratio;
 xnew =  xprime';
end

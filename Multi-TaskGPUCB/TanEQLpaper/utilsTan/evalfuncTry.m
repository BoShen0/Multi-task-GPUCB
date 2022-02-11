function Y=evalfunc(D)
%Shubert Function + Schwefel Function (*)
u1=1*(D(:,1)+D(:,3))+1; v1=1*(D(:,2)+D(:,4))+1;
n=size(u1,1); 
u2=50*(D(:,1)+D(:,3))+50; v2=50*(D(:,2)+D(:,4))+50;
Y=[(sum(repmat(1:5,n,1).*cos(repmat((1:5)+1,n,1).*repmat(u1,1,5)+repmat((1:5),n,1)),2).*sum(repmat(1:5,n,1).*cos(repmat((1:5)+1,n,1).*repmat(v1,1,5)+repmat((1:5),n,1)),2))'];
%     (418.9829*2-u2.*sin(sqrt(abs(u2)))-v2.*sin(sqrt(abs(v2))))']./[3.72858741286651;57.5662387791625];
end

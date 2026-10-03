function [x,w]=gauss01(n)
% n-point Gauss-Legendre rule on [0,1].
validateattributes(n,{'numeric'},{'scalar','integer','positive'});
if n==1, x=.5; w=1; return; end
k=(1:n-1)'; b=k./sqrt(4*k.^2-1);
[V,D]=eig(diag(b,1)+diag(b,-1));
[x,idx]=sort(diag(D)); w=(V(1,idx).^2)'; x=(x+1)/2;
end

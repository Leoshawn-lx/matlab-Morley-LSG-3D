function [b,w]=triangle_quadrature(order)
% Normalized triangle weights. order=2 uses the exact degree-two midpoint rule.
if order==2,b=[0 .5 .5;.5 0 .5;.5 .5 0];w=ones(3,1)/3;return;end
[q,g]=gauss01(order); [U,V]=ndgrid(q,q); [A,B]=ndgrid(g,g);
u=U(:);v=V(:);b=[1-u,u.*(1-v),u.*v];w=2*A(:).*B(:).*u;
end

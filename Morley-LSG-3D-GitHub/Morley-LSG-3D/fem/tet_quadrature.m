function [b,w]=tet_quadrature(order)
% Weights normalized to sum one; multiply by physical tetrahedron volume.
if order==2
    a=.5854101966249685; c=.1381966011250105;
    b=c*ones(4,4)+(a-c)*eye(4); w=ones(4,1)/4;return;
end
[q,g]=gauss01(order); [U,V,W]=ndgrid(q,q,q); [A,B,C]=ndgrid(g,g,g);
u=U(:);v=V(:);t=W(:);
b=[1-u,u.*(1-v),u.*v.*(1-t),u.*v.*t];
w=6*A(:).*B(:).*C(:).*u.^2.*v;
end

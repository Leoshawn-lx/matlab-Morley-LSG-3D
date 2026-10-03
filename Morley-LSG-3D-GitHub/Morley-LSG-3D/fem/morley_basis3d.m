function [phi,grad,H]=morley_basis3d(xyz,g)
% phi(n,10), grad(n,10,3); H(10,3,3), constant on an affine tetrahedron.
q=bsxfun(@rdivide,bsxfun(@minus,xyz,g.center),g.scale);
if nargout==1, phi=p2_monomials(q)*g.C; return; end
[m,d,hm]=p2_monomials(q); phi=m*g.C; n=size(xyz,1); grad=zeros(n,10,3);
for j=1:3,grad(:,:,j)=(d(:,:,j)*g.C)/g.scale;end
if nargout<3,return;end
H=zeros(10,3,3);
for i=1:3
    for j=1:3,H(:,i,j)=(hm(:,i,j)'*g.C)'/g.scale^2;end
end
end

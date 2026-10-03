function [m,d,H]=p2_monomials(q)
% [1,x,y,z,x^2,y^2,z^2,xy,xz,yz], derivatives with respect to q.
x=q(:,1); y=q(:,2); z=q(:,3); n=size(q,1);
m=[ones(n,1),x,y,z,x.^2,y.^2,z.^2,x.*y,x.*z,y.*z];
if nargout<2,return;end
d=zeros(n,10,3); d(:,2,1)=1; d(:,3,2)=1; d(:,4,3)=1;
d(:,5,1)=2*x; d(:,6,2)=2*y; d(:,7,3)=2*z;
d(:,8,1)=y; d(:,8,2)=x; d(:,9,1)=z; d(:,9,3)=x; d(:,10,2)=z; d(:,10,3)=y;
if nargout<3,return;end
H=zeros(10,3,3); H(5,1,1)=2; H(6,2,2)=2; H(7,3,3)=2;
H(8,1,2)=1; H(8,2,1)=1; H(9,1,3)=1; H(9,3,1)=1; H(10,2,3)=1; H(10,3,2)=1;
end

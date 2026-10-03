function g=tet_geometry(v,normals)
% Six edge averages, four normal-derivative averages in GLOBAL normals.
% The face rows are scaled by h before inversion to avoid unit imbalance.
e=[1 2;1 3;1 4;2 3;2 4;3 4]; f=[2 3 4;1 4 3;1 2 4;1 3 2];
g.center=mean(v,1); g.scale=max(sqrt(sum((v(e(:,1),:)-v(e(:,2),:)).^2,2)));
detB=det(v(2:4,:)-v(1,:));
if ~isfinite(detB) || abs(detB)<1e-13*g.scale^3,error('annular:degenerateTet','Degenerate or extremely flat tetrahedron.');end
g.volume=abs(detB)/6;
if nargin<2
    normals=zeros(4,3);
    for k=1:4
        p=v(f(k,:),:); nn=cross(p(2,:)-p(1,:),p(3,:)-p(1,:)); nn=nn/norm(nn);
        if dot(nn,v(k,:)-mean(p,1))>0,nn=-nn;end
        normals(k,:)=nn;
    end
end
q=bsxfun(@rdivide,bsxfun(@minus,v,g.center),g.scale); D=zeros(10,10);
for k=1:6
    p=[q(e(k,1),:);mean(q(e(k,:),:),1);q(e(k,2),:)];
    D(k,:)=[1 4 1]*p2_monomials(p)/6;
end
for k=1:4
    [~,dm]=p2_monomials(mean(q(f(k,:),:),1));
    D(k+6,:)=(reshape(dm,10,3)*normals(k,:)')';
end
g.localRcond=rcond(D);
if g.localRcond<1e-13,error('annular:singularBasis','Local Morley moment matrix is nearly singular.');end
g.C=D\diag([ones(1,6),g.scale*ones(1,4)]);
end

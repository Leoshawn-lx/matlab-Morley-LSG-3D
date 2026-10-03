function [J,S,rd,fd,brow,Tmap]=local_face_js(m,f,lambda,mu)
%LOCAL_FACE_JS Projected JS face operator.
% J = J(v,u)-J(u,v), with rows=test DOFs and columns=trial DOFs.
% Tmap maps a constant boundary displacement vector gbar to the local
% vector J(v,gbar)=|F| T' gbar.  It is used only for the nonhomogeneous
% boundary-data RHS correction on boundary faces.
cells=double(m.face2elem(f,:));cells=cells(cells>0);ns=numel(cells);
n=m.faceNormal(f,:)';[b,w]=triangle_quadrature(2);x=b*m.node(double(m.face(f,:)),:);
U=zeros(3,30*ns);T=U;d=zeros(1,30*ns);rd=zeros(1,30*ns);fd=rd;
for k=1:ns
    t=cells(k);g=element_geometry(m,t);pm=w'*morley_basis3d(x,g);
    [~,dg]=morley_basis3d(g.center,g);G=reshape(dg,10,3);
    sg=1;if k==2,sg=-1;end
    [rd((k-1)*30+(1:30)),fd((k-1)*30+(1:30))]=element_dofs(m,t);
    for c=1:3
        idx=(k-1)*30+(c-1)*10+(1:10);ec=zeros(3,1);ec(c)=1;
        U(c,idx)=sg*pm;
        T(:,idx)=(mu*(ec*(G*n)'+n(c)*G')+lambda*n*G(:,c)')/ns;
        d(idx)=sg*G(:,c)';
    end
end
Tmap=m.faceArea(f)*T';
J=Tmap*U-m.faceArea(f)*(U'*T);
S=zeros(size(J));brow=zeros(size(d));
if ns==2
    brow=lambda*sqrt(m.faceH(f)*m.faceArea(f))*d;
    S=brow'*brow;
end
end

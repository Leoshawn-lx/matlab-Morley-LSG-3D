function m=build_morley3d_topology(r,cfg)
% Tagged boundary triangles define constrained EDGES and FACES, not vertices.
node=double(r.node); elem=double(r.tet); nt=size(elem,1);
e=[1 2;1 3;1 4;2 3;2 4;3 4]; f=[2 3 4;1 4 3;1 2 4;1 3 2];
req={'Inner','Outer','Bottom','Top','Domain'}; tags=zeros(1,5);
for k=1:5
    q=find(strcmp({r.physicalNames.name},req{k}));
    if numel(q)~=1,error('annular:physicalGroups','Missing or duplicate group %s.',req{k});end
    tags(k)=r.physicalNames(q).tag;
    if r.physicalNames(q).dim~=2+(k==5),error('annular:physicalGroups','Wrong group dimension.');end
end
if isempty(r.tri) || any(r.tetTag~=tags(5)),error('annular:physicalGroups','Missing triangles or invalid Domain tags.');end
for k=1:nt
    v=node(elem(k,:),:); d=det(v(2:4,:)-v(1,:));
    h=max(sqrt(sum((v(e(:,1),:)-v(e(:,2),:)).^2,2)));
    if ~isfinite(d) || abs(d)<1e-13*h^3,error('annular:degenerateTet','Invalid tetrahedron %d.',k);end
    if d<0,elem(k,[2 3])=elem(k,[3 2]);end
end
te=zeros(6*nt,2); tf=zeros(4*nt,3);
for k=1:6,te((k-1)*nt+(1:nt),:)=sort(elem(:,e(k,:)),2);end
for k=1:4,tf((k-1)*nt+(1:nt),:)=sort(elem(:,f(k,:)),2);end
[edge,~,ie]=unique(te,'rows'); [face,~,iff]=unique(tf,'rows');
e2e=reshape(ie,nt,6); e2f=reshape(iff,nt,4); nf=size(face,1); ne=size(edge,1);
counts=accumarray(iff,1,[nf,1]);
if any(counts>2),error('annular:nonManifold','A face has more than two incident tetrahedra.');end
[~,ord]=sort(iff); first=[1;1+cumsum(counts(1:end-1))]; owner=repmat((1:nt)',4,1);
f2e=zeros(nf,2); f2e(:,1)=owner(ord(first)); inside=counts==2;
f2e(inside,2)=owner(ord(first(inside)+1)); bnd=find(~inside);
[ok,tagFace]=ismember(sort(r.tri,2),face,'rows');
if ~all(ok) || numel(unique(tagFace))~=numel(tagFace) || ~isequal(sort(tagFace),bnd)
    error('annular:boundaryTags','Tagged triangles must cover every topological boundary face exactly once.');
end
bt=zeros(nf,1); phys=zeros(nf,1); phys(tagFace)=r.triTag;
for k=1:4
    ids=tagFace(r.triTag==tags(k));
    if isempty(ids),error('annular:physicalGroups','Empty surface group %s.',req{k});end
    bt(ids)=k;
end
if any(bt(bnd)==0),error('annular:boundaryTags','Unknown boundary physical tag.');end
be=unique(sort([face(bnd,[1 2]);face(bnd,[1 3]);face(bnd,[2 3])],2),'rows');
[ok,beid]=ismember(be,edge,'rows'); assert(all(ok));
scalarNdof=ne+nf; fixedScalar=unique([beid;ne+bnd]);
freeScalar=setdiff((1:scalarNdof)',fixedScalar); nfs=numel(freeScalar);
free=[freeScalar;freeScalar+scalarNdof;freeScalar+2*scalarNdof];
if numel(free)>cfg.maxFreeDof,error('annular:dofGuard','%d free DOFs exceed maxFreeDof=%d. Change config deliberately to proceed.',numel(free),cfg.maxFreeDof);end
if numel(free)>cfg.warnFreeDof,warning('annular:largeMesh','Large mesh: %d free DOFs; factorization memory may be substantial.',numel(free));end
ctr=(node(elem(:,1),:)+node(elem(:,2),:)+node(elem(:,3),:)+node(elem(:,4),:))/4;
p1=node(face(:,1),:);p2=node(face(:,2),:);p3=node(face(:,3),:);
nn=cross(p2-p1,p3-p1,2); area=sqrt(sum(nn.^2,2))/2;
if any(area<=0),error('annular:degenerateFace','Zero-area face.');end
nn=bsxfun(@rdivide,nn,2*area);
flip=sum(nn.*(ctr(f2e(:,1),:)-(p1+p2+p3)/3),2)>0; nn(flip,:)=-nn(flip,:);
hface=max([sqrt(sum((p2-p1).^2,2)),sqrt(sum((p3-p1).^2,2)),sqrt(sum((p3-p2).^2,2))],[],2);
C=zeros(10,10,nt); vol=zeros(nt,1); h=zeros(nt,1); rc=zeros(nt,1);
for k=1:nt
    g=tet_geometry(node(elem(k,:),:),nn(e2f(k,:),:));
    C(:,:,k)=g.C; vol(k)=g.volume; h(k)=g.scale; rc(k)=g.localRcond;
end
m.node=node; m.elem=uint32(elem); m.edge=uint32(edge);m.face=uint32(face);
m.elem2edge=uint32(e2e);m.elem2face=uint32(e2f);m.face2elem=uint32(f2e);
m.faceNormal=nn;m.faceArea=area;m.faceH=hface;m.faceCount=counts;
m.boundaryFaces=uint32(bnd);m.boundaryEdges=uint32(beid);m.boundaryType=uint8(bt);m.physicalTag=phys;
m.scalarNdof=scalarNdof;m.vectorNdof=3*scalarNdof;m.nFree=numel(free);m.nFreeScalar=nfs;
m.freeDof=free;m.fixedDof=[fixedScalar;fixedScalar+scalarNdof;fixedScalar+2*scalarNdof];
m.scalarMap=uint32([e2e,ne+e2f]); m.scalarFreeMap=zeros(scalarNdof,1,'uint32');
m.scalarFreeMap(freeScalar)=uint32(1:nfs);
m.center=ctr;m.scale=h;m.volume=vol;m.C=C;m.localRcond=rc;m.hmax=max(h);
m.totalVolume=sum(vol);m.signature=cfg_signature(cfg);
[~,uid]=fileparts(tempname);m.uid=uid;
end

function [gFull,diag] = build_pullout_boundary_dofs(m,ref,cfg)
%BUILD_PULLOUT_BOUNDARY_DOFS Interpolate the nonhomogeneous essential data
%into the 3-D generalized-Morley boundary moments.
%
% Scalar Morley DOFs:
%   1:ne       edge averages of u;
%   ne+1:...   face averages of d_n u in the GLOBAL face normal stored by m.
%
% The exact physical boundary is circular while the mesh is chordal.
% Displacement values on Inner/Outer are imposed from the physical boundary
% constants ub/0.  On Top/Bottom, w(r) is averaged along boundary edges.
% For a cylindrical boundary face, grad u at the radial projection is
% dotted with the planar mesh normal n_h before face averaging.  Hence the
% data converge to the exact d_n trace as the polygonal cylinder is refined.

ne=size(m.edge,1); nf=size(m.face,1); ns=m.scalarNdof;
gScalar=zeros(ns,1); assigned=false(ns,1);

bnd=double(m.boundaryFaces(:)); bt=double(m.boundaryType(:));
% Edge-to-boundary-group incidence.
F=double(m.face(bnd,:));
allE=sort([F(:,[1 2]);F(:,[1 3]);F(:,[2 3])],2);
[ok,eid]=ismember(allE,double(m.edge),'rows');
if ~all(ok), error('stage4:edgeMap','A boundary-face edge is absent from the global edge table.'); end
etype=[bt(bnd);bt(bnd);bt(bnd)];
edgeGroup=false(ne,4);
for k=1:numel(eid), edgeGroup(eid(k),etype(k))=true; end
bedge=double(m.boundaryEdges(:));
if any(~any(edgeGroup(bedge,:),2)), error('stage4:edgeGroup','Unclassified boundary edge.'); end
if any(edgeGroup(bedge,1) & edgeGroup(bedge,2)), error('stage4:edgeGroup','An edge cannot be both Inner and Outer.'); end

[q1,w1]=gauss01(cfg.edgeQuadOrder);
for kk=1:numel(bedge)
    e=bedge(kk); p0=m.node(double(m.edge(e,1)),:); p1=m.node(double(m.edge(e,2)),:);
    if edgeGroup(e,1)
        val=cfg.ub;
    elseif edgeGroup(e,2)
        val=0;
    else
        x=(1-q1).*p0 + q1.*p1;
        r=hypot(x(:,1),x(:,2));
        [ww,~,~,~]=eval_aifantis_reference_1d(ref,r);
        val=sum(w1.*ww);
    end
    gScalar(e)=val; assigned(e)=true;
end

% Face-normal derivative moments.
[bq,wq]=triangle_quadrature(cfg.faceQuadOrder);
faceVal=zeros(numel(bnd),1);
for kk=1:numel(bnd)
    f=bnd(kk); typ=bt(f); sid=ne+f;
    if typ==3 || typ==4
        val=0; % z=0,H; reference is z-independent.
    else
        x=bq*m.node(double(m.face(f,:)),:);
        rr=hypot(x(:,1),x(:,2));
        er=[x(:,1)./rr,x(:,2)./rr,zeros(size(rr))];
        nh=m.faceNormal(f,:);
        if typ==1, wp=ref.wp(1); else, wp=ref.wp(end); end
        dnu=wp*(er*nh(:));
        val=sum(wq.*dnu);
    end
    gScalar(sid)=val; assigned(sid)=true; faceVal(kk)=val;
end

fixedScalar=unique([bedge; ne+bnd]);
if ~all(assigned(fixedScalar)), error('stage4:boundaryInterpolation','Not every fixed scalar DOF received data.'); end
if any(abs(gScalar(setdiff((1:ns)',fixedScalar)))>0), error('stage4:boundaryInterpolation','Interior scalar DOF was assigned boundary data.'); end

% u1=u2=0; only u3 is nonzero.
gFull=zeros(m.vectorNdof,1);
gFull(2*ns+(1:ns))=gScalar;
if any(abs(gFull(double(m.freeDof)))>0), error('stage4:boundaryInterpolation','A free DOF received prescribed data.'); end

% Diagnostics by physical group.
diag=struct();
diag.nBoundaryEdges=numel(bedge); diag.nBoundaryFaces=numel(bnd);
diag.nNonzeroFixed=nnz(gFull(double(m.fixedDof)));
diag.maxAbsBoundaryValue=max(abs(gFull(double(m.fixedDof))));
diag.innerEdgeRange=range_or_zero(gScalar(bedge(edgeGroup(bedge,1))));
diag.outerEdgeRange=range_or_zero(gScalar(bedge(edgeGroup(bedge,2))));
diag.planeOnlyEdgeRange=range_or_zero(gScalar(bedge(~edgeGroup(bedge,1)&~edgeGroup(bedge,2))));
diag.faceMean=zeros(1,4);diag.faceMin=zeros(1,4);diag.faceMax=zeros(1,4);diag.nFace=zeros(1,4);
for t=1:4
    v=faceVal(bt(bnd)==t);diag.nFace(t)=numel(v);
    if ~isempty(v),diag.faceMean(t)=mean(v);diag.faceMin(t)=min(v);diag.faceMax(t)=max(v);end
end
diag.innerExactNormalTrace=ref.trace.inner.dnu3;
diag.outerExactNormalTrace=ref.trace.outer.dnu3;
diag.planeExactNormalTrace=0;
diag.innerFaceMeanRelDiff=abs(diag.faceMean(1)-diag.innerExactNormalTrace)/max(abs(diag.innerExactNormalTrace),eps);
diag.outerFaceMeanRelDiff=abs(diag.faceMean(2)-diag.outerExactNormalTrace)/max(abs(diag.outerExactNormalTrace),eps);
end

function a=range_or_zero(v)
if isempty(v),a=[0,0];else,a=[min(v),max(v)];end
end

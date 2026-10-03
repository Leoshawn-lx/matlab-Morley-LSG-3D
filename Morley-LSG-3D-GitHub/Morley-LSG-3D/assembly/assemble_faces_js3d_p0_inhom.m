function [A,b,stats]=assemble_faces_js3d_p0_inhom(m,lambda,gFull,ref,cfg,A,b)
%ASSEMBLE_FACES_JS3D_P0_INHOM Add projected JS face terms and nonzero-BC RHS.
%
% The manuscript scheme is derived for homogeneous clamped data.  On a
% boundary face, the nonsymmetric consistency pair is the analogue of
%   -<sigma(u)n,v> + <sigma(v)n,u>.
% For prescribed displacement g, consistency therefore requires the RHS
% term +<sigma(v)n,g>, i.e. +J_boundary(v,g).  Merely eliminating the
% prescribed Morley DOFs from the homogeneous matrix misses this term.
%
% We keep the bilinear form A unchanged and add only the boundary-data RHS
% correction.  Interior J and S are exactly unchanged.
if nargin<7,b=zeros(m.nFree,1);end
if nargin<6,A=sparse(m.nFree,m.nFree);end
nf=size(m.face,1);chunk=cfg.chunkFaces;n=m.nFree;timer=tic;
stats.nJFaces=nf;stats.nSFaces=nnz(m.face2elem(:,2)>0);
stats.maxJSkewDefect=0;stats.maxSSymmetryDefect=0;
bJg=zeros(n,1);
for start=1:chunk:nf
    stop=min(nf,start+chunk-1);cap=3600*(stop-start+1);
    I=zeros(cap,1);JJ=I;V=I;pos=0;
    for f=start:stop
        [J,S,rd,fd,~,Tmap]=local_face_js(m,f,lambda,cfg.mu); K=J+S;
        stats.maxJSkewDefect=max(stats.maxJSkewDefect,norm(J+J','fro')/max(norm(J,'fro'),realmin));
        stats.maxSSymmetryDefect=max(stats.maxSSymmetryDefect,norm(S-S','fro')/max(norm(S,'fro'),realmin));
        gloc=gFull(fd(:));
        bloc=-K*gloc;

        % Nonhomogeneous boundary consistency:
        % RHS += J_boundary(v,g) = |F| (P0 sigma(v)n) . mean_F(g).
        if m.face2elem(f,2)==0
            gbar=boundary_displacement_mean(m,f,ref,cfg);
            dataLocal=Tmap*gbar;
            bloc=bloc+dataLocal;
            fr=find(rd>0);
            if ~isempty(fr), bJg=accum_rhs(bJg,rd(fr),dataLocal(fr)); end
        end

        fr=find(rd>0); if ~isempty(fr), b=accum_rhs(b,rd(fr),bloc(fr)); end
        [ii,jj]=ndgrid(rd,rd);keep=ii>0 & jj>0 & K~=0;count=nnz(keep);idx=pos+(1:count);
        I(idx)=ii(keep);JJ(idx)=jj(keep);V(idx)=K(keep);pos=pos+count;
    end
    A=A+sparse(I(1:pos),JJ(1:pos),V(1:pos),n,n);
    if start==1 || stop==nf || mod(start-1,chunk*100)==0,fprintf('  faces %d/%d\n',stop,nf);end
end
stats.seconds=toc(timer);
stats.boundaryDataRhsNorm=norm(bJg);
end

function gbar=boundary_displacement_mean(m,f,ref,cfg)
typ=double(m.boundaryType(f));
if typ==1
    gbar=[0;0;cfg.ub];
elseif typ==2
    gbar=[0;0;0];
elseif typ==3 || typ==4
    [bq,wq]=triangle_quadrature(cfg.faceQuadOrder);
    x=bq*m.node(double(m.face(f,:)),:);
    r=hypot(x(:,1),x(:,2));
    [ww,~,~,~]=eval_aifantis_reference_1d(ref,r);
    gbar=[0;0;sum(wq.*ww)];
else
    error('stage4b:boundaryType','Boundary face %d has invalid type %d.',f,typ);
end
end

function b=accum_rhs(b,idx,val)
b=b+accumarray(idx(:),val(:),[numel(b),1],@sum,0);
end

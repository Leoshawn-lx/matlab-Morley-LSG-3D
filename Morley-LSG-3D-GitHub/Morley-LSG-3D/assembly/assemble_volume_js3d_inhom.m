function [A,b,stats]=assemble_volume_js3d_inhom(m,lambda,iota,gFull,cfg)
% Reduced free/free volume matrix plus lifting RHS -K_FB*g_B.
nt=size(m.elem,1); n=m.nFree; A=sparse(n,n); b=zeros(n,1);
chunk=cfg.chunkElements; timer=tic;
for start=1:chunk:nt
    stop=min(nt,start+chunk-1); cap=900*(stop-start+1);
    I=zeros(cap,1);JJ=I;V=I;pos=0;
    for k=start:stop
        [rd,fd]=element_dofs(m,k);
        K=local_volume_js(m.node(double(m.elem(k,:)),:),element_geometry(m,k),lambda,iota,cfg.mu);
        gloc=gFull(fd(:)); bloc=-K*gloc;
        fr=find(rd>0); if ~isempty(fr), b=accum_rhs(b,rd(fr),bloc(fr)); end
        [ii,jj]=ndgrid(rd,rd);keep=ii>0 & jj>0 & K~=0;count=nnz(keep);idx=pos+(1:count);
        I(idx)=ii(keep);JJ(idx)=jj(keep);V(idx)=K(keep);pos=pos+count;
    end
    A=A+sparse(I(1:pos),JJ(1:pos),V(1:pos),n,n);
    if start==1 || stop==nt || mod(start-1,chunk*20)==0,fprintf('  volume %d/%d\n',stop,nt);end
end
stats.seconds=toc(timer);stats.symmetryDefect=norm(A-A','fro')/max(norm(A,'fro'),realmin);
end

function b=accum_rhs(b,idx,val)
b=b+accumarray(idx(:),val(:),[numel(b),1],@sum,0);
end

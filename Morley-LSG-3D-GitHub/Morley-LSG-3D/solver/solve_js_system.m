function s=solve_js_system(A,b,cfg,nFreeScalar,x0)
%SOLVE_JS_SYSTEM Scaled nonsymmetric JS solve with component-block preconditioning.
%
% The reduced vector unknowns are stored componentwise as
%   [u1_free ; u2_free ; u3_free].
% For N=8 the full symmetric-part ICT can break down numerically even though
% the JS energy is coercive.  Here we avoid factorizing the full 79056-by-79056
% symmetric part.  Instead, for each displacement component c=1,2,3 we form
% the principal block
%   H_c = 0.5*(A_cc + A_cc'),
% apply a symmetric AMD ordering, and build an incomplete Cholesky factor.
% If ICT still breaks down for a block, that block automatically falls back
% to a symmetric Gauss-Seidel (SGS) triangular preconditioner.
%
% The finite element matrix A is NOT altered.  All shifts/approximations are
% used only inside the preconditioner.  Acceptance is always based on
%   ||A*x-b||/||b|| <= cfg.solverTol.

timer=tic; n=size(A,1);
if size(A,2)~=n || numel(b)~=n
    error('pullout:systemShape','Invalid linear-system shape.');
end
if nargin<4 || isempty(nFreeScalar)
    if mod(n,3)~=0
        error('pullout:blockSize','Need nFreeScalar, or n divisible by 3.');
    end
    nFreeScalar=n/3;
end
if n~=3*nFreeScalar
    error('pullout:blockSize','Expected n = 3*nFreeScalar.');
end
if any(~isfinite(nonzeros(A))) || any(~isfinite(b))
    error('pullout:systemData','Nonfinite system entries.');
end
if nargin<5 || isempty(x0)
    x0=zeros(n,1);
end
if numel(x0)~=n || any(~isfinite(x0))
    error('pullout:initialGuess','Initial guess must be a finite vector of system size.');
end

s=struct('x',zeros(n,1),'flag',1,'rawFlag',1,'relres',Inf,'accepted',false, ...
    'solver','','iterations',[],'seconds',0,'factorizationSeconds',0, ...
    'preconditionerSeconds',0,'krylovSeconds',0,'preconditionedResidual',Inf, ...
    'preconditionerNnz',0,'blockMode',strings(3,1),'blockNnz',zeros(3,1), ...
    'blockDiagComp',nan(3,1),'blockDropTol',nan(3,1), ...
    'refinementResidualHistory',[],'refinementCorrectionRel',[], ...
    'refinementFlags',[],'refinementIterations',[]);
if norm(b)==0
    s.flag=0; s.rawFlag=0; s.relres=0; s.accepted=true;
    s.solver='zero-rhs'; s.seconds=toc(timer); return;
end

% Symmetric diagonal scaling.  This leaves the linear system equivalent but
% makes every diagonal entry approximately one, which helps the Krylov solve.
ad=full(diag(A));
if any(~isfinite(ad)) || any(ad<=0)
    error('pullout:diagonal','Expected a positive finite diagonal before scaling.');
end
d=1./sqrt(ad);
As=spdiags(d,0,n,n)*A*spdiags(d,0,n,n);
bs=d.*b;

if strcmp(cfg.solver,'direct')
    fprintf('  Building sparse LU for systemNdof=%d ...\n',n);
    tFac=tic;
    if exist('decomposition','class') || exist('decomposition','file')
        fac=decomposition(As,'lu'); apply=@(rhs)fac\rhs;
    else
        [L,U,P,Q,R]=lu(As); apply=@(rhs)Q*(U\(L\(P*(R\rhs))));
    end
    s.factorizationSeconds=toc(tFac);
    fprintf('  Sparse LU ready in %.2f s.\n',s.factorizationSeconds);
    y=apply(bs); x=d.*y;
    s.rawFlag=0; s.solver='scaled-direct-LU'; s.iterations=0;
else
    if ~strcmp(cfg.solver,'gmres-block')
        error('pullout:solver','Stage 4D expects cfg.solver=''gmres-block'' or ''direct''.');
    end

    fprintf('  Building 3 component-block preconditioners from symmetric principal blocks ...\n');
    tPre=tic;
    blocks=repmat(struct('mode','','L',[],'p',[],'T',[],'diagv',[]),3,1);
    totalPnnz=0;
    for c=1:3
        idx=(c-1)*nFreeScalar+(1:nFreeScalar);
        Ac=As(idx,idx);
        Hc=0.5*(Ac+Ac');
        Hc=0.5*(Hc+Hc');
        clear Ac
        symDef=norm(Hc-Hc','fro')/max(norm(Hc,'fro'),eps);
        fprintf('    block %d: size=%d, nnz(Hc)=%d, symmetry defect %.3e\n', ...
            c,nFreeScalar,nnz(Hc),symDef);
        dc0=full(diag(Hc));
        if any(~isfinite(dc0)) || any(dc0<=0)
            error('pullout:blockDiagonal','Block %d has nonpositive/nonfinite diagonal.',c);
        end

        % Reordering is essential: the previous full-matrix natural-order ICT
        % broke down with a negative pivot.
        p=symamd(Hc);
        Hp=Hc(p,p);
        clear Hc

        built=false;
        dropList=cfg.blockIcholDropTolList;
        compList=cfg.blockIcholDiagCompList;
        micholList={'off','on'};
        for idr=1:numel(dropList)
            for im=1:numel(micholList)
                for idc=1:numel(compList)
                    opts=struct('type','ict','droptol',dropList(idr), ...
                        'diagcomp',compList(idc),'michol',micholList{im});
                    fprintf('      ICT try: dropTol=%.1e, diagComp=%.1e, michol=%s ...\n', ...
                        dropList(idr),compList(idc),micholList{im});
                    try
                        L=ichol(Hp,opts);
                        built=true;
                        blocks(c).mode='ICT';
                        blocks(c).L=L;
                        blocks(c).p=p(:);
                        s.blockDiagComp(c)=compList(idc);
                        s.blockDropTol(c)=dropList(idr);
                        s.blockNnz(c)=nnz(L);
                        totalPnnz=totalPnnz+nnz(L);
                        fprintf('      block %d ICT ready: nnz(L)=%d\n',c,nnz(L));
                        break
                    catch ME
                        fprintf('        failed: %s\n',ME.message);
                    end
                end
                if built, break; end
            end
            if built, break; end
        end

        if ~built
            % Guaranteed low-setup-cost fallback.  Hp has positive diagonal;
            % T=tril(Hp) is therefore nonsingular.  M=T*D^{-1}*T' is used as
            % a symmetric Gauss-Seidel preconditioner.  No factorization is
            % performed, so this fallback cannot suffer an IC negative pivot.
            fprintf('      ICT failed for block %d; using SGS triangular fallback.\n',c);
            T=tril(Hp);
            blocks(c).mode='SGS';
            blocks(c).T=T;
            blocks(c).diagv=full(diag(Hp));
            blocks(c).p=p(:);
            s.blockNnz(c)=nnz(T);
            totalPnnz=totalPnnz+nnz(T);
        end
        s.blockMode(c)=string(blocks(c).mode);
        clear Hp p L T
    end
    s.preconditionerNnz=totalPnnz;
    s.preconditionerSeconds=toc(tPre);
    fprintf('  Block preconditioner ready in %.2f s; total stored nnz=%d\n', ...
        s.preconditionerSeconds,totalPnnz);
    fprintf('  block modes: %s / %s / %s\n', ...
        s.blockMode(1),s.blockMode(2),s.blockMode(3));

    Mfun=@(r)apply_component_blocks(r,blocks,nFreeScalar);
    fprintf('  Starting GMRES (restart=%d, tol=%.1e, maxit=%d) ...\n', ...
        cfg.gmresRestart,cfg.gmresTol,cfg.gmresMaxit);
    tKrylov=tic;
    y0=x0./d;
    initialTrueResidual=norm(A*x0-b)/norm(b);
    fprintf('  Initial-guess true relative residual = %.3e\n',initialTrueResidual);
    [y,flag,preRes,it]=gmres(As,bs,cfg.gmresRestart,cfg.gmresTol, ...
        cfg.gmresMaxit,Mfun,[],y0);
    s.rawFlag=flag;
    s.preconditionedResidual=preRes;
    s.iterations=it;
    s.solver='scaled-GMRES-component-block-refined';
    fprintf('  Initial GMRES finished: flag=%d; reported relres=%.3e; iter=[%d %d]\n', ...
        flag,preRes,it(1),it(2));

    % Residual-correction iterative refinement.  This is important for the
    % lambda=1e5 system because ||b|| is O(1e10): a small relative residual
    % can still leave an O(1) absolute residual in weakly controlled modes.
    x=d.*y;
    hist=zeros(cfg.refinementSteps+1,1);
    corr=zeros(cfg.refinementSteps,1);
    rflags=zeros(cfg.refinementSteps,1);
    rits=zeros(cfg.refinementSteps,2);
    hist(1)=norm(A*x-b)/norm(b);
    fprintf('  True residual before refinement = %.3e\n',hist(1));
    for ir=1:cfg.refinementSteps
        if hist(ir)<=cfg.solverTol
            hist=hist(1:ir);
            corr=corr(1:ir-1);
            rflags=rflags(1:ir-1);
            rits=rits(1:ir-1,:);
            break
        end
        r=b-A*x;
        rs=d.*r;
        fprintf('  Refinement %d: solving correction, ||r||/||b||=%.3e ...\n',ir,hist(ir));
        [dy,fr,rr,itr]=gmres(As,rs,cfg.gmresRestart,cfg.refinementTol, ...
            cfg.refinementMaxit,Mfun,[],zeros(n,1));
        dx=d.*dy;
        corr(ir)=norm(dx)/max(norm(x),eps);
        y=y+dy;
        x=d.*y;
        hist(ir+1)=norm(A*x-b)/norm(b);
        rflags(ir)=fr;
        rits(ir,:)=itr;
        fprintf('    correction flag=%d, reported relres=%.3e, iter=[%d %d], rel ||dx||=%.3e\n', ...
            fr,rr,itr(1),itr(2),corr(ir));
        fprintf('    true residual after refinement %d = %.3e\n',ir,hist(ir+1));
    end
    s.refinementResidualHistory=hist(:).';
    s.refinementCorrectionRel=corr(:).';
    s.refinementFlags=rflags(:).';
    s.refinementIterations=rits;
    s.krylovSeconds=toc(tKrylov);
    fprintf('  Total Krylov+refinement time %.2f s\n',s.krylovSeconds);
end

s.x=x;
s.relres=norm(A*x-b)/norm(b);
s.seconds=toc(timer);
s.accepted=all(isfinite(x)) && isfinite(s.relres) && s.relres<=cfg.solverTol;
s.flag=double(~s.accepted);
fprintf('  True unpreconditioned relative residual = %.3e\n',s.relres);
end

function z=apply_component_blocks(r,blocks,nb)
% Return M^{-1}r for a block-diagonal component preconditioner.
z=zeros(size(r));
for c=1:3
    idx=(c-1)*nb+(1:nb);
    rc=r(idx);
    p=blocks(c).p;
    rp=rc(p);
    if strcmp(blocks(c).mode,'ICT')
        L=blocks(c).L;
        zp=L'\(L\rp);
    else
        T=blocks(c).T;
        dv=blocks(c).diagv;
        zp=T'\(dv.*(T\rp));
    end
    zc=zeros(nb,1);
    zc(p)=zp;
    z(idx)=zc;
end
end

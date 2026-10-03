function s=solve_fullH_gmres(A,b,cfg,x0,mode)
%SOLVE_FULLH_GMRES GMRES with full coupled symmetric-part preconditioner.
% mode = 'chol' : exact sparse Cholesky of H=(As+As')/2 (N=4 calibration only)
% mode = 'ict'  : incomplete Cholesky of the same full coupled H
%
% IMPORTANT: A and b are never changed.  Any shift is used only in the
% preconditioner.  The returned true residual is measured on the original
% unscaled system.

n=size(A,1);
if nargin<4 || isempty(x0), x0=zeros(n,1); end
if nargin<5, mode='ict'; end
if norm(b)==0
    s=struct('x',zeros(n,1),'relres',0,'accepted',true,'flag',0,'iter',[0 0], ...
        'reportedRelres',0,'setupSeconds',0,'solveSeconds',0,'mode',mode, ...
        'shift',0,'nnzFactor',0,'scaledResidual',0);
    return
end

ad=full(diag(A));
if any(~isfinite(ad)) || any(ad<=0)
    error('stage5d:diag','Expected positive finite diagonal for symmetric scaling.');
end
d=1./sqrt(ad);
As=spdiags(d,0,n,n)*A*spdiags(d,0,n,n);
bs=d.*b;
y0=x0./d;
H=0.5*(As+As');
H=0.5*(H+H');

fprintf('  Full coupled H: size=%d, nnz=%d, symmetry defect %.3e\n', ...
    n,nnz(H),norm(H-H','fro')/max(norm(H,'fro'),eps));
p=symamd(H);
Hp=H(p,p);
clear H

setupT=tic;
shiftUsed=0;
if strcmpi(mode,'chol')
    fprintf('  Building exact sparse Cholesky of full H (N=4 calibration) ...\n');
    [L,pflag]=chol(Hp,'lower');
    if pflag~=0
        base=max(mean(abs(full(diag(Hp)))),1);
        shiftList=[1e-14,1e-12,1e-10,1e-8,1e-6];
        built=false;
        for k=1:numel(shiftList)
            tau=shiftList(k)*base;
            fprintf('    chol retry with preconditioner-only shift %.3e ...\n',tau);
            [L,pflag]=chol(Hp+tau*speye(n),'lower');
            if pflag==0
                shiftUsed=tau; built=true; break
            end
        end
        if ~built
            error('stage5d:chol','Full-H Cholesky failed even after tiny shifts.');
        end
    end
elseif strcmpi(mode,'ict')
    fprintf('  Building full coupled ICT preconditioner ...\n');
    built=false;
    dropList=cfg.fullHDropTolList;
    compList=cfg.fullHDiagCompList;
    micholList={'on','off'};
    for im=1:numel(micholList)
        for idr=1:numel(dropList)
            for idc=1:numel(compList)
                opts=struct('type','ict','droptol',dropList(idr), ...
                    'diagcomp',compList(idc),'michol',micholList{im});
                fprintf('    ICT try: dropTol=%.1e, diagComp=%.1e, michol=%s ...\n', ...
                    dropList(idr),compList(idc),micholList{im});
                try
                    L=ichol(Hp,opts);
                    shiftUsed=compList(idc);
                    built=true;
                    fprintf('    full-H ICT ready: nnz(L)=%d\n',nnz(L));
                    break
                catch ME
                    fprintf('      failed: %s\n',ME.message);
                end
            end
            if built, break; end
        end
        if built, break; end
    end
    if ~built
        error('stage5d:ict','Could not build a full-H ICT preconditioner.');
    end
else
    error('stage5d:mode','Unknown mode %s.',mode);
end
setupSeconds=toc(setupT);

Mfun=@(r)apply_factor(r,L,p);
restart=cfg.fullHGmresRestart;
tol=cfg.fullHGmresTol;
maxit=cfg.fullHGmresMaxit;
initialOrig=norm(A*x0-b)/norm(b);
initialScaled=norm(As*y0-bs)/norm(bs);
fprintf('  GMRES(%s): initial original/scaled residual %.3e / %.3e\n', ...
    mode,initialOrig,initialScaled);
solveT=tic;
[y,flag,relres,it]=gmres(As,bs,restart,tol,maxit,Mfun,[],y0);
solveSeconds=toc(solveT);
x=d.*y;
origRes=norm(A*x-b)/norm(b);
scaledRes=norm(As*y-bs)/norm(bs);
fprintf('  GMRES(%s) done: flag=%d, reported=%.3e, iter=[%d %d]\n', ...
    mode,flag,relres,it(1),it(2));
fprintf('  true original/scaled residual = %.3e / %.3e\n',origRes,scaledRes);

s=struct();
s.x=x;
s.relres=origRes;
s.scaledResidual=scaledRes;
s.accepted=all(isfinite(x)) && origRes<=cfg.fullHAcceptTol;
s.flag=flag;
s.iter=it;
s.reportedRelres=relres;
s.setupSeconds=setupSeconds;
s.solveSeconds=solveSeconds;
s.mode=mode;
s.shift=shiftUsed;
s.nnzFactor=nnz(L);
end

function z=apply_factor(r,L,p)
rp=r(p);
zp=L'\(L\rp);
z=zeros(size(r));
z(p)=zp;
end

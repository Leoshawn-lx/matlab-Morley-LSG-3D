function s=solve_energy_transform_gmres(A,b,cfg,x0,label)
%SOLVE_ENERGY_TRANSFORM_GMRES Two-sided exact-H transformed GMRES.
%
% The FE system is NOT changed.  We first use symmetric diagonal scaling
%   As = D*A*D,  bs = D*b,  x = D*y,
% and form the exact symmetric part
%   H = (As+As')/2.
% With Hp = H(p,p) = L*L' after a symmetric AMD permutation, set
%   y_p = L'\z.
% The transformed system is
%   B*z = c,
%   B = L\(Ap/L'),  c = L\bp,
% whose symmetric part is the identity in exact arithmetic.  This keeps the
% full u1-u2-u3 lambda coupling in the preconditioner/coordinate transform.
%
% A memory guard uses symbfact before the numerical Cholesky factorization.

if nargin<4 || isempty(x0), x0=zeros(size(A,1),1); end
if nargin<5, label=''; end
n=size(A,1);
if norm(b)==0
    s=struct('x',zeros(n,1),'relres',0,'scaledResidual',0,'transformedResidual',0, ...
        'accepted',true,'flag',0,'iter',[0 0],'reportedRelres',0,'setupSeconds',0, ...
        'solveSeconds',0,'nnzFactor',0,'estimatedNnzFactor',0,'estimatedFactorGB',0, ...
        'skipped',false,'skipReason','');
    return
end

ad=full(diag(A));
if any(~isfinite(ad)) || any(ad<=0)
    error('stage5e:diag','Expected positive finite diagonal for symmetric scaling.');
end
d=1./sqrt(ad);
D=spdiags(d,0,n,n);
As=D*A*D;
bs=d.*b;
y0=x0./d;
clear D

H=0.5*(As+As');
H=0.5*(H+H');
fprintf('  [%s] scaled H: size=%d, nnz=%d, symmetry defect %.3e\n', ...
    label,n,nnz(H),norm(H-H','fro')/max(norm(H,'fro'),eps));

p=symamd(H);
Hp=H(p,p);
Ap=As(p,p);
bp=bs(p);
y0p=y0(p);
clear H As bs y0

% Symbolic fill estimate before numerical Cholesky.
setupT=tic;
fprintf('  [%s] estimating exact-H Cholesky fill with symbfact ...\n',label);
try
    counts=symbfact(Hp);
    estNnz=sum(double(counts));
catch ME
    warning('stage5e:symbfact','symbfact failed (%s); proceeding without fill estimate.',ME.message);
    estNnz=NaN;
end
if isfinite(estNnz)
    % Conservative rough storage estimate: values + row indices + workspace allowance.
    estGB=24*estNnz/2^30;
    fprintf('  [%s] estimated nnz(L)=%.3e, rough factor storage %.2f GB\n',label,estNnz,estGB);
else
    estGB=NaN;
end

availGB=NaN;
if ispc
    try
        [~,sys]=memory;
        availGB=double(sys.PhysicalMemory.Available)/2^30;
        fprintf('  [%s] Windows-reported available physical memory %.2f GB\n',label,availGB);
    catch
    end
end
limitGB=cfg.energyMaxEstimatedFactorGB;
if isfinite(availGB)
    limitGB=min(limitGB,cfg.energyMaxAvailFraction*availGB);
end
if isfinite(estGB) && estGB>limitGB
    reason=sprintf('estimated factor storage %.2f GB exceeds guard %.2f GB',estGB,limitGB);
    fprintf('  [%s] MEMORY GUARD: %s. Cholesky is NOT attempted.\n',label,reason);
    s=struct('x',[],'relres',Inf,'scaledResidual',Inf,'transformedResidual',Inf, ...
        'accepted',false,'flag',99,'iter',[0 0],'reportedRelres',Inf, ...
        'setupSeconds',toc(setupT),'solveSeconds',0,'nnzFactor',0, ...
        'estimatedNnzFactor',estNnz,'estimatedFactorGB',estGB,'skipped',true, ...
        'skipReason',reason);
    return
end

fprintf('  [%s] building exact sparse Cholesky of full coupled H ...\n',label);
[L,pflag]=chol(Hp,'lower');
shiftUsed=0;
if pflag~=0
    base=max(mean(abs(full(diag(Hp)))),1);
    shiftList=[1e-14,1e-13,1e-12,1e-11,1e-10,1e-9,1e-8];
    built=false;
    for k=1:numel(shiftList)
        tau=shiftList(k)*base;
        fprintf('    [%s] chol retry with preconditioner-only shift %.3e ...\n',label,tau);
        [L,pflag]=chol(Hp+tau*speye(n),'lower');
        if pflag==0
            shiftUsed=tau; built=true; break
        end
    end
    if ~built
        error('stage5e:chol','Exact full-H Cholesky failed even after tiny shifts.');
    end
end
setupSeconds=toc(setupT);
fprintf('  [%s] Cholesky ready: nnz(L)=%d, setup %.2f s, shift %.3e\n', ...
    label,nnz(L),setupSeconds,shiftUsed);
clear Hp

% Two-sided energy-coordinate transform.
Bfun=@(z) L\(Ap*(L'\z));
c=L\bp;
z0=L'*y0p;
initialOrig=norm(A*x0-b)/norm(b);
initialTrans=norm(Bfun(z0)-c)/norm(c);
fprintf('  [%s] initial original/transformed residual %.3e / %.3e\n', ...
    label,initialOrig,initialTrans);

solveT=tic;
[z,flag,reported,it]=gmres(Bfun,c,cfg.energyRestart,cfg.energyTol, ...
    cfg.energyMaxit,[],[],z0);
solveSeconds=toc(solveT);
yp=L'\z;
y=zeros(n,1); y(p)=yp;
x=d.*y;
origRes=norm(A*x-b)/norm(b);
% Re-evaluate scaled residual without rebuilding the full scaled matrix.
scaledR=d.*(A*x-b);
scaledB=d.*b;
scaledRes=norm(scaledR)/norm(scaledB);
transRes=norm(Bfun(z)-c)/norm(c);
fprintf('  [%s] transformed GMRES: flag=%d, reported=%.3e, iter=[%d %d], time %.2f s\n', ...
    label,flag,reported,it(1),it(2),solveSeconds);
fprintf('  [%s] true original/scaled/transformed residual = %.3e / %.3e / %.3e\n', ...
    label,origRes,scaledRes,transRes);

s=struct();
s.x=x;
s.relres=origRes;
s.scaledResidual=scaledRes;
s.transformedResidual=transRes;
s.accepted=all(isfinite(x)) && transRes<=cfg.energyAcceptTransformedTol && origRes<=cfg.energyAcceptOriginalTol;
s.flag=flag;
s.iter=it;
s.reportedRelres=reported;
s.setupSeconds=setupSeconds;
s.solveSeconds=solveSeconds;
s.nnzFactor=nnz(L);
s.estimatedNnzFactor=estNnz;
s.estimatedFactorGB=estGB;
s.shift=shiftUsed;
s.skipped=false;
s.skipReason='';
end

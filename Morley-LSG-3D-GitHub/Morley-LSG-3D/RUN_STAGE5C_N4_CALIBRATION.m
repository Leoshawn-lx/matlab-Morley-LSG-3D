clear; clc;
cfg=setup_stage5();
iota=cfg.iota;
lambdaLo=1;
lambdaHi=1e5;

fprintf('==============================================================\n');
fprintf('Stage 5C-v3: FAIR N=4 solver calibration for lambda=1e5\n');
fprintf('Finite element formulation is UNCHANGED.\n');
fprintf('Purpose: compare Stage-5B-style GMRES against direct LU on N=4.\n');
fprintf('IMPORTANT: the lambda=1 solution is reused as the lambda=1e5 GMRES initial guess,\n');
fprintf('exactly matching the strategy used in Stage 5B on N=8.\n');
fprintf('mu=%g, iota=%g, lambdaLo=%g, lambdaHi=%g\n',cfg.mu,iota,lambdaLo,lambdaHi);
fprintf('==============================================================\n\n');

pref=stage5_reference_config(cfg);
ref=solve_aifantis_reference_1d(pref,iota);
fprintf('1-D reference max BC residual = %.3e\n',ref.maxBoundaryResidual);
fprintf('inner/outer d_n u3 = %.9g / %.9g\n\n',ref.trace.inner.dnu3,ref.trace.outer.dnu3);

N=4;
m=build_stage4_mesh(N,cfg);
[g,bd]=build_pullout_boundary_dofs(m,ref,cfg);
fprintf('Mesh N=4: tet=%d, DOF=%d, free=%d, hmax=%.9g\n', ...
    size(m.elem,1),m.vectorNdof,m.nFree,m.hmax);

%% Step 1: lambda=1 direct solve, used ONLY as the Stage-5B-style initial guess.
fprintf('\n--------------------------------------------------------------\n');
fprintf('Step 1: N=4 lambda=1 direct solve for initial guess\n');
fprintf('--------------------------------------------------------------\n');
[A1,b1,vs1]=assemble_volume_js3d_inhom(m,lambdaLo,iota,g,cfg);
[A1,b1,fs1]=assemble_faces_js3d_p0_inhom(m,lambdaLo,g,ref,cfg,A1,b1);
cfgD=cfg; cfgD.solver='direct';
sol1=solve_js_system(A1,b1,cfgD,m.nFreeScalar,zeros(m.nFree,1));
if ~sol1.accepted, error('stage5c:lambda1direct','N=4 lambda=1 direct solve failed.'); end
x0=sol1.x;
fprintf('lambda=1 direct true residual = %.3e\n',sol1.relres);
clear A1 b1 vs1 fs1

%% Step 2: lambda=1e5 direct LU reference.
fprintf('\n--------------------------------------------------------------\n');
fprintf('Step 2: N=4 lambda=1e5 DIRECT reference solve\n');
fprintf('--------------------------------------------------------------\n');
[A,b,vs]=assemble_volume_js3d_inhom(m,lambdaHi,iota,g,cfg);
[A,b,fs]=assemble_faces_js3d_p0_inhom(m,lambdaHi,g,ref,cfg,A,b);
fprintf('System: size=%d, nnz=%d, ||b||=%.6e\n',size(A,1),nnz(A),norm(b));
fprintf('Checks: volume sym %.3e, J skew %.3e, S sym %.3e\n', ...
    vs.symmetryDefect,fs.maxJSkewDefect,fs.maxSSymmetryDefect);
fprintf('Boundary J(v,g) RHS norm = %.6e\n',fs.boundaryDataRhsNorm);

solD=solve_js_system(A,b,cfgD,m.nFreeScalar,zeros(m.nFree,1));
if ~solD.accepted, error('stage5c:direct','N=4 lambda=1e5 direct solve failed.'); end
uD=g; uD(double(m.freeDof))=solD.x;
errD=compute_pullout_errors(m,uD,ref,iota,cfg);
pathD=compute_radial_profile(m,uD,ref,cfg);
fprintf('DIRECT errors: L2=%.6e, strain=%.6e, gradStrain=%.6e, threebar=%.6e\n', ...
    errD.L2,errD.strain,errD.gradStrain,errD.threebar);
fprintf('               H1=%.6e, H2=%.6e, ||div uh||=%.6e\n', ...
    errD.H1,errD.H2,errD.divUh);
fprintf('DIRECT profile max error = %.6e\n',pathD.maxAbsError);

%% Step 3: Stage-5B-equivalent GMRES solve, same low-lambda initial guess.
fprintf('\n--------------------------------------------------------------\n');
fprintf('Step 3: N=4 lambda=1e5 Stage-5B-style GMRES calibration\n');
fprintf('--------------------------------------------------------------\n');
cfgG=cfg;
cfgG.solver='gmres-block';
% Match the Stage-5B policy rather than the over-strict previous Stage-5C run.
cfgG.solverTol=1e-8;
cfgG.gmresTol=1e-9;
cfgG.gmresRestart=40;
cfgG.gmresMaxit=200;
cfgG.refinementSteps=0;
fprintf('Using converged lambda=1 solution as GMRES initial guess.\n');
solG=solve_js_system(A,b,cfgG,m.nFreeScalar,x0);

% IMPORTANT: this is a diagnostic.  Even if the iterative residual does not
% meet 1e-8, continue and compare the actual FE vector/errors with direct LU.
if ~solG.accepted
    warning('stage5c:diagnosticOnly', ...
        'GMRES did not meet %.1e (actual %.3e). Continuing comparison for diagnosis only.', ...
        cfgG.solverTol,solG.relres);
end
uG=g; uG(double(m.freeDof))=solG.x;
errG=compute_pullout_errors(m,uG,ref,iota,cfg);
pathG=compute_radial_profile(m,uG,ref,cfg);

coefDiff=norm(solG.x-solD.x)/max(norm(solD.x),eps);
fullDiff=norm(uG-uD)/max(norm(uD),eps);
compDiff=zeros(3,1);
for c=1:3
    idx=(c-1)*m.nFreeScalar+(1:m.nFreeScalar);
    compDiff(c)=norm(solG.x(idx)-solD.x(idx))/max(norm(solD.x(idx)),eps);
end

fprintf('\n==============================================================\n');
fprintf('N=4 solver calibration result\n');
fprintf('==============================================================\n');
fprintf('Direct residual = %.6e\n',solD.relres);
fprintf('GMRES residual  = %.6e (accepted by 1e-8 criterion: %d)\n',solG.relres,solG.accepted);
fprintf('GMRES-vs-direct free-vector relative difference = %.6e\n',coefDiff);
fprintf('GMRES-vs-direct full-vector relative difference = %.6e\n',fullDiff);
fprintf('Component relative differences u1/u2/u3 = %.6e / %.6e / %.6e\n', ...
    compDiff(1),compDiff(2),compDiff(3));
fprintf('GMRES errors: L2=%.6e, strain=%.6e, gradStrain=%.6e, threebar=%.6e\n', ...
    errG.L2,errG.strain,errG.gradStrain,errG.threebar);
fprintf('              H1=%.6e, H2=%.6e, ||div uh||=%.6e\n', ...
    errG.H1,errG.H2,errG.divUh);
fprintf('Error ratios GMRES/direct:\n');
fprintf('  L2=%.8f, strain=%.8f, gradStrain=%.8f, threebar=%.8f\n', ...
    errG.L2/errD.L2,errG.strain/errD.strain, ...
    errG.gradStrain/errD.gradStrain,errG.threebar/errD.threebar);
fprintf('  H1=%.8f, H2=%.8f, div ratio=%.6e\n', ...
    errG.H1/errD.H1,errG.H2/errD.H2,errG.divUh/max(errD.divUh,eps));
fprintf('Profile max errors direct / GMRES = %.6e / %.6e\n', ...
    pathD.maxAbsError,pathG.maxAbsError);
fprintf('Max radial-profile difference GMRES vs direct = %.6e\n', ...
    max(abs(pathG.uh3-pathD.uh3)));

% Interpretation flags.  Do not call a failed residual solve "accepted";
% these flags only classify whether the FE quantities are sensitive to the
% iterative error at N=4.
errRatio=[errG.L2/errD.L2,errG.strain/errD.strain, ...
          errG.gradStrain/errD.gradStrain,errG.threebar/errD.threebar];
maxErrRatioDev=max(abs(errRatio-1));
if maxErrRatioDev<=1e-3 && coefDiff<=1e-5
    fprintf('DIAGNOSTIC VERDICT: GMRES and direct LU agree very closely on N=4.\n');
    fprintf('This argues AGAINST solver error as the main source of the Stage-5B N=8 high-order anomaly.\n');
elseif maxErrRatioDev<=1e-2 && coefDiff<=1e-4
    fprintf('DIAGNOSTIC VERDICT: GMRES and direct LU agree reasonably on N=4.\n');
    fprintf('Solver contamination appears small but is not completely excluded.\n');
else
    fprintf('DIAGNOSTIC VERDICT: GMRES differs materially from direct LU on N=4.\n');
    fprintf('Do NOT trust the N=8 high-order lambda comparison until the large-lambda solver is improved.\n');
end
fprintf('==============================================================\n');

result=struct();
result.cfg=cfg;
result.direct=struct('solver',solD,'error',errD,'path',pathD);
result.gmres=struct('solver',solG,'error',errG,'path',pathG);
result.coefDiff=coefDiff;
result.fullDiff=fullDiff;
result.compDiff=compDiff;
result.errorRatios=errRatio;
result.maxErrorRatioDeviation=maxErrRatioDev;
save(fullfile(cfg.resultDir,'stage5c_v3_N4_fair_calibration.mat'),'result','-v7.3');
fprintf('Saved diagnostic MAT: %s\n',fullfile(cfg.resultDir,'stage5c_v3_N4_fair_calibration.mat'));

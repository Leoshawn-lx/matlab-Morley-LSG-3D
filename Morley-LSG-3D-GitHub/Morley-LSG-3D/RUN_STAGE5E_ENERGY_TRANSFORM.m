clear; clc;
cfg=setup_stage5();
iota=cfg.iota;
lambdaLo=1;
lambdaHi=1e5;

% Exact full-H two-sided energy transform settings.
cfg.energyRestart=80;
cfg.energyTol=1e-10;
cfg.energyMaxit=120;
cfg.energyAcceptTransformedTol=5e-9;
cfg.energyAcceptOriginalTol=1e-8;
% Guard against an accidental out-of-memory Cholesky on N=8.
cfg.energyMaxEstimatedFactorGB=6.0;
cfg.energyMaxAvailFraction=0.45;

fprintf('==============================================================\n');
fprintf('Stage 5E: energy-coordinate exact-full-H solver validation\n');
fprintf('Finite element formulation is UNCHANGED.\n');
fprintf('Part A validates the solver on N=4 against direct LU.\n');
fprintf('Part B runs N=8 lambda=1e5 only if Part A passes and the\n');
fprintf('symbolic Cholesky memory guard is satisfied.\n');
fprintf('mu=%g, iota=%g, lambdaLo=%g, lambdaHi=%g\n',cfg.mu,iota,lambdaLo,lambdaHi);
fprintf('==============================================================\n\n');

pref=stage5_reference_config(cfg);
ref=solve_aifantis_reference_1d(pref,iota);
fprintf('1-D reference max BC residual = %.3e\n',ref.maxBoundaryResidual);
fprintf('inner/outer d_n u3 = %.9g / %.9g\n',ref.trace.inner.dnu3,ref.trace.outer.dnu3);

%% =============================================================
% Part A: N=4 direct calibration
%% =============================================================
fprintf('\n==============================================================\n');
fprintf('PART A: N=4 direct-LU calibration of energy-transform GMRES\n');
fprintf('==============================================================\n');
N=4;
m4=build_stage4_mesh(N,cfg);
[g4,~]=build_pullout_boundary_dofs(m4,ref,cfg);

% lambda=1 direct initial guess + errors
fprintf('\n--- N=4 lambda=1 direct ---\n');
[A4lo,b4lo]=assemble_volume_js3d_inhom(m4,lambdaLo,iota,g4,cfg);
[A4lo,b4lo]=assemble_faces_js3d_p0_inhom(m4,lambdaLo,g4,ref,cfg,A4lo,b4lo);
cfgD=cfg; cfgD.solver='direct'; cfgD.solverTol=1e-10;
sol4lo=solve_js_system(A4lo,b4lo,cfgD,m4.nFreeScalar,zeros(m4.nFree,1));
if ~sol4lo.accepted, error('stage5e:n4lo','N=4 lambda=1 direct solve failed.'); end
u4lo=g4; u4lo(double(m4.freeDof))=sol4lo.x;
err4lo=compute_pullout_errors(m4,u4lo,ref,iota,cfg);
clear A4lo b4lo u4lo

% lambda=1e5 direct reference
fprintf('\n--- N=4 lambda=1e5 direct reference ---\n');
[A4,b4]=assemble_volume_js3d_inhom(m4,lambdaHi,iota,g4,cfg);
[A4,b4]=assemble_faces_js3d_p0_inhom(m4,lambdaHi,g4,ref,cfg,A4,b4);
sol4D=solve_js_system(A4,b4,cfgD,m4.nFreeScalar,zeros(m4.nFree,1));
if ~sol4D.accepted, error('stage5e:n4direct','N=4 lambda=1e5 direct solve failed.'); end
u4D=g4; u4D(double(m4.freeDof))=sol4D.x;
err4D=compute_pullout_errors(m4,u4D,ref,iota,cfg);
fprintf('N=4 DIRECT high-lambda: L2=%.6e, strain=%.6e, gradStrain=%.6e, threebar=%.6e, H2=%.6e, div=%.6e\n', ...
    err4D.L2,err4D.strain,err4D.gradStrain,err4D.threebar,err4D.H2,err4D.divUh);

% energy-transform exact H
fprintf('\n--- N=4 lambda=1e5 exact-H energy transform ---\n');
sol4E=solve_energy_transform_gmres(A4,b4,cfg,sol4lo.x,'N4-lam1e5');
if sol4E.skipped, error('stage5e:n4skip','Unexpected N=4 memory-guard skip: %s',sol4E.skipReason); end
u4E=g4; u4E(double(m4.freeDof))=sol4E.x;
err4E=compute_pullout_errors(m4,u4E,ref,iota,cfg);
coef4=norm(sol4E.x-sol4D.x)/max(norm(sol4D.x),eps);
ratio4=[err4E.L2/err4D.L2,err4E.strain/err4D.strain,err4E.gradStrain/err4D.gradStrain,err4E.threebar/err4D.threebar];
fprintf('\nN=4 ENERGY-TRANSFORM CALIBRATION\n');
fprintf('  free-vector relative diff = %.6e\n',coef4);
fprintf('  error ratios energy/direct = [%.8f %.8f %.8f %.8f]\n',ratio4);
fprintf('  H2 ratio=%.8f, div direct/energy = %.6e / %.6e\n',err4E.H2/err4D.H2,err4D.divUh,err4E.divUh);

% Primary criterion is reproduction of FE error functionals, not a needlessly
% strict coefficient-vector tolerance in an extremely ill-conditioned system.
pass4=max(abs(ratio4-1))<=5e-3 && abs(err4E.H2/err4D.H2-1)<=5e-3 && coef4<=2e-3;
if pass4
    fprintf('  PART A VERDICT: PASS. Exact-H energy transform reproduces N=4 direct-LU FE errors.\n');
else
    fprintf('  PART A VERDICT: FAIL. Do not run N=8.\n');
    save(fullfile(cfg.resultDir,'stage5e_partA_failed.mat'),'sol4D','sol4E','err4D','err4E','coef4','ratio4','cfg','-v7.3');
    return
end
clear A4 b4 u4D u4E

%% =============================================================
% Part B: N=8 lambda robustness with trusted lambda=1 + energy high-lambda
%% =============================================================
fprintf('\n==============================================================\n');
fprintf('PART B: N=8 lambda=1 vs lambda=1e5\n');
fprintf('==============================================================\n');
N=8;
m8=build_stage4_mesh(N,cfg);
[g8,bd8]=build_pullout_boundary_dofs(m8,ref,cfg);
fprintf('Mesh N=8: tet=%d, DOF=%d, free=%d, hmax=%.9g\n',size(m8.elem,1),m8.vectorNdof,m8.nFree,m8.hmax);

% lambda=1: validated component-block solver used only to provide the baseline
% FE solution and high-lambda initial guess.
fprintf('\n--- N=8 lambda=1 validated component-block solve ---\n');
[A8lo,b8lo]=assemble_volume_js3d_inhom(m8,lambdaLo,iota,g8,cfg);
[A8lo,b8lo]=assemble_faces_js3d_p0_inhom(m8,lambdaLo,g8,ref,cfg,A8lo,b8lo);
cfgLo=cfg; cfgLo.solver='gmres-block'; cfgLo.solverTol=1e-8; cfgLo.gmresTol=1e-9; ...
    cfgLo.gmresRestart=40; cfgLo.gmresMaxit=200; cfgLo.refinementSteps=0;
sol8lo=solve_js_system(A8lo,b8lo,cfgLo,m8.nFreeScalar,zeros(m8.nFree,1));
if ~sol8lo.accepted, error('stage5e:n8lo','N=8 lambda=1 baseline solve failed.'); end
u8lo=g8; u8lo(double(m8.freeDof))=sol8lo.x;
err8lo=compute_pullout_errors(m8,u8lo,ref,iota,cfg);
path8lo=compute_radial_profile(m8,u8lo,ref,cfg);
fprintf('N=8 lambda=1: L2=%.6e, strain=%.6e, gradStrain=%.6e, threebar=%.6e, H2=%.6e, div=%.6e\n', ...
    err8lo.L2,err8lo.strain,err8lo.gradStrain,err8lo.threebar,err8lo.H2,err8lo.divUh);
clear A8lo b8lo u8lo

% lambda=1e5 exact-H energy transformed solve.
fprintf('\n--- N=8 lambda=1e5 exact-H energy-transform solve ---\n');
[A8,b8,vs8]=assemble_volume_js3d_inhom(m8,lambdaHi,iota,g8,cfg);
[A8,b8,fs8]=assemble_faces_js3d_p0_inhom(m8,lambdaHi,g8,ref,cfg,A8,b8);
fprintf('N=8 high-lambda system: size=%d, nnz=%d, ||b||=%.6e\n',size(A8,1),nnz(A8),norm(b8));
fprintf('Checks: volume sym %.3e, J skew %.3e, S sym %.3e\n',vs8.symmetryDefect,fs8.maxJSkewDefect,fs8.maxSSymmetryDefect);
sol8hi=solve_energy_transform_gmres(A8,b8,cfg,sol8lo.x,'N8-lam1e5');
if sol8hi.skipped
    fprintf('\nPART B STOPPED SAFELY BY MEMORY GUARD: %s\n',sol8hi.skipReason);
    save(fullfile(cfg.resultDir,'stage5e_memory_guard.mat'),'sol8hi','cfg','-v7.3');
    return
end
if ~sol8hi.accepted
    warning('stage5e:n8accept','N=8 energy-transform solver did not meet residual acceptance; results are diagnostic only.');
end
u8hi=g8; u8hi(double(m8.freeDof))=sol8hi.x;
err8hi=compute_pullout_errors(m8,u8hi,ref,iota,cfg);
path8hi=compute_radial_profile(m8,u8hi,ref,cfg);
clear A8 b8 u8hi

ratio8=[err8hi.L2/err8lo.L2,err8hi.strain/err8lo.strain,err8hi.gradStrain/err8lo.gradStrain,err8hi.threebar/err8lo.threebar];
profileDiff=max(abs(path8hi.uh3-path8lo.uh3));
profileDiffRel=profileDiff/max(max(abs(ref.w)),cfg.ub);
dh=log(m4.hmax/m8.hmax);
rateLo=[log(err4lo.L2/err8lo.L2),log(err4lo.strain/err8lo.strain),log(err4lo.gradStrain/err8lo.gradStrain),log(err4lo.threebar/err8lo.threebar)]/dh;
rateHi=[log(err4D.L2/err8hi.L2),log(err4D.strain/err8hi.strain),log(err4D.gradStrain/err8hi.gradStrain),log(err4D.threebar/err8hi.threebar)]/dh;

fprintf('\n==============================================================\n');
fprintf('Stage 5E FINAL N=8 lambda diagnostic\n');
fprintf('==============================================================\n');
fprintf('lambda=1 errors    : L2=%.6e, strain=%.6e, gradStrain=%.6e, threebar=%.6e\n', ...
    err8lo.L2,err8lo.strain,err8lo.gradStrain,err8lo.threebar);
fprintf('lambda=1e5 errors  : L2=%.6e, strain=%.6e, gradStrain=%.6e, threebar=%.6e\n', ...
    err8hi.L2,err8hi.strain,err8hi.gradStrain,err8hi.threebar);
fprintf('high/low ratios    : [%.6f %.6f %.6f %.6f]\n',ratio8);
fprintf('N=4->8 rates lam=1 : [%.3f %.3f %.3f %.3f]\n',rateLo);
fprintf('N=4->8 rates 1e5   : [%.3f %.3f %.3f %.3f]\n',rateHi);
fprintf('divUh lam=1 / 1e5  : %.6e / %.6e\n',err8lo.divUh,err8hi.divUh);
fprintf('radial profile diff relative to max reference displacement = %.3e\n',profileDiffRel);
fprintf('energy solver residuals original/scaled/transformed = %.3e / %.3e / %.3e\n', ...
    sol8hi.relres,sol8hi.scaledResidual,sol8hi.transformedResidual);
fprintf('exact-H factor nnz=%d (symbolic estimate %.3e, rough %.2f GB)\n', ...
    sol8hi.nnzFactor,sol8hi.estimatedNnzFactor,sol8hi.estimatedFactorGB);
fprintf('==============================================================\n');

summary=table([lambdaLo;lambdaHi],[err8lo.L2;err8hi.L2],[err8lo.strain;err8hi.strain], ...
    [err8lo.gradStrain;err8hi.gradStrain],[err8lo.threebar;err8hi.threebar], ...
    [err8lo.H1;err8hi.H1],[err8lo.H2;err8hi.H2],[err8lo.divUh;err8hi.divUh], ...
    [rateLo(1);rateHi(1)],[rateLo(2);rateHi(2)],[rateLo(3);rateHi(3)],[rateLo(4);rateHi(4)], ...
    'VariableNames',{'lambda','L2','strain','gradStrain','threebar','H1','H2','divUh', ...
    'rateL2','rateStrain','rateGradStrain','rateThreebar'});
writetable(summary,fullfile(cfg.resultDir,'stage5e_N8_energy_transform.csv'));
save(fullfile(cfg.resultDir,'stage5e_N8_energy_transform.mat'),'summary','err4lo','err4D','err4E', ...
    'err8lo','err8hi','sol8hi','ratio4','ratio8','rateLo','rateHi','profileDiffRel','bd8','cfg','-v7.3');

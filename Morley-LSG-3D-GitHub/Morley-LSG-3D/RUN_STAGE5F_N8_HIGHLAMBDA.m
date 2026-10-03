clear; clc;
cfg=setup_stage5();
iota=cfg.iota;
lambdaHi=1e5;

% Exact full-H energy transform settings validated on N=4 in Stage 5E.
cfg.energyRestart=80;
cfg.energyTol=1e-10;
cfg.energyMaxit=120;
cfg.energyAcceptTransformedTol=5e-9;
cfg.energyAcceptOriginalTol=1e-8;
cfg.energyMaxEstimatedFactorGB=6.0;
cfg.energyMaxAvailFraction=0.45;

fprintf('==============================================================\n');
fprintf('Stage 5F-B: FRESH-SESSION N=8 lambda=1e5 exact-H solve\n');
fprintf('Finite element formulation is UNCHANGED.\n');
fprintf('Run this script only in a FRESH MATLAB session.\n');
fprintf('The N=8 lambda=1 initial guess must first be prepared by\n');
fprintf('RUN_STAGE5F_PREPARE_X0.\n');
fprintf('==============================================================\n\n');

% Early memory check before spending time on assembly.
if ispc
    try
        [~,sys]=memory;
        avail0=double(sys.PhysicalMemory.Available)/2^30;
        fprintf('Windows-reported available physical memory at start: %.2f GB\n',avail0);
        if avail0 < 7.5
            fprintf('\nMEMORY PRECHECK STOP: less than 7.5 GB is currently available.\n');
            fprintf('Close memory-heavy applications, restart MATLAB, and run ONLY this script again.\n');
            return
        end
    catch
    end
end

x0File=fullfile(cfg.resultDir,'stage5f_N8_lambda1_x0.mat');
if ~exist(x0File,'file')
    error('stage5f:x0missing','Missing %s. Run RUN_STAGE5F_PREPARE_X0 first.',x0File);
end
S=load(x0File,'x0','err8lo','path8lo','hmax8','nFree8');
x0=S.x0;
err8lo=S.err8lo;
path8lo=S.path8lo;
hmax8=S.hmax8;
nFree8=S.nFree8;
clear S

% Verified N=4 direct-LU baselines from Stage 5A.
baseFile=fullfile(cfg.resultDir,'stage5a_verified_N4_baseline.csv');
T4=readtable(baseFile);
rowHi=find(T4.lambda==lambdaHi,1);
rowLo=find(T4.lambda==1,1);
if isempty(rowHi) || isempty(rowLo)
    error('stage5f:baseline','Verified N=4 baseline CSV is incomplete.');
end

pref=stage5_reference_config(cfg);
ref=solve_aifantis_reference_1d(pref,iota);
N=8;
m8=build_stage4_mesh(N,cfg);
[g8,bd8]=build_pullout_boundary_dofs(m8,ref,cfg);
if m8.nFree~=nFree8 || abs(m8.hmax-hmax8)>1e-12
    error('stage5f:mesh','Prepared lambda=1 initial guess does not match the rebuilt N=8 mesh.');
end

fprintf('Mesh N=8: tet=%d, DOF=%d, free=%d, hmax=%.9g\n', ...
    size(m8.elem,1),m8.vectorNdof,m8.nFree,m8.hmax);
fprintf('\n--- N=8 lambda=1e5 exact-H energy-transform solve ---\n');
[A,b,vs]=assemble_volume_js3d_inhom(m8,lambdaHi,iota,g8,cfg);
[A,b,fs]=assemble_faces_js3d_p0_inhom(m8,lambdaHi,g8,ref,cfg,A,b);
fprintf('System: size=%d, nnz=%d, ||b||=%.6e\n',size(A,1),nnz(A),norm(b));
fprintf('Checks: volume sym %.3e, J skew %.3e, S sym %.3e\n', ...
    vs.symmetryDefect,fs.maxJSkewDefect,fs.maxSSymmetryDefect);
clear vs fs

sol=solve_energy_transform_gmres(A,b,cfg,x0,'N8-lam1e5-fresh');
if sol.skipped
    fprintf('\nStage 5F-B STOPPED SAFELY BY MEMORY GUARD: %s\n',sol.skipReason);
    fprintf('No numerical high-lambda conclusion should be drawn from this run.\n');
    save(fullfile(cfg.resultDir,'stage5f_memory_guard.mat'),'sol','cfg','-v7.3');
    return
end
if ~sol.accepted
    warning('stage5f:accept','Energy-transform solve did not satisfy acceptance thresholds; results are diagnostic only.');
end

u=g8;
u(double(m8.freeDof))=sol.x;
err8hi=compute_pullout_errors(m8,u,ref,iota,cfg);
path8hi=compute_radial_profile(m8,u,ref,cfg);
clear A b u

ratio8=[err8hi.L2/err8lo.L2,err8hi.strain/err8lo.strain, ...
    err8hi.gradStrain/err8lo.gradStrain,err8hi.threebar/err8lo.threebar];
profileDiff=max(abs(path8hi.uh3-path8lo.uh3));
profileDiffRel=profileDiff/max(max(abs(ref.w)),cfg.ub);
dh=log(T4.hmax(rowLo)/m8.hmax);
rateLo=[log(T4.L2(rowLo)/err8lo.L2),log(T4.strain(rowLo)/err8lo.strain), ...
    log(T4.gradStrain(rowLo)/err8lo.gradStrain),log(T4.threebar(rowLo)/err8lo.threebar)]/dh;
rateHi=[log(T4.L2(rowHi)/err8hi.L2),log(T4.strain(rowHi)/err8hi.strain), ...
    log(T4.gradStrain(rowHi)/err8hi.gradStrain),log(T4.threebar(rowHi)/err8hi.threebar)]/dh;

fprintf('\n==============================================================\n');
fprintf('Stage 5F FINAL N=8 lambda diagnostic\n');
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
    sol.relres,sol.scaledResidual,sol.transformedResidual);
fprintf('exact-H factor nnz=%d (symbolic estimate %.3e, rough %.2f GB)\n', ...
    sol.nnzFactor,sol.estimatedNnzFactor,sol.estimatedFactorGB);
fprintf('==============================================================\n');

summary=table([1;lambdaHi],[err8lo.L2;err8hi.L2],[err8lo.strain;err8hi.strain], ...
    [err8lo.gradStrain;err8hi.gradStrain],[err8lo.threebar;err8hi.threebar], ...
    [err8lo.H1;err8hi.H1],[err8lo.H2;err8hi.H2],[err8lo.divUh;err8hi.divUh], ...
    [rateLo(1);rateHi(1)],[rateLo(2);rateHi(2)],[rateLo(3);rateHi(3)],[rateLo(4);rateHi(4)], ...
    'VariableNames',{'lambda','L2','strain','gradStrain','threebar','H1','H2','divUh', ...
    'rateL2','rateStrain','rateGradStrain','rateThreebar'});
writetable(summary,fullfile(cfg.resultDir,'stage5f_N8_exactH_fresh.csv'));
save(fullfile(cfg.resultDir,'stage5f_N8_exactH_fresh.mat'),'summary','err8lo','err8hi','sol', ...
    'ratio8','rateLo','rateHi','profileDiffRel','bd8','cfg','-v7.3');

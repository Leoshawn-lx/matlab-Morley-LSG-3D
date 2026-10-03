clear; clc;
cfg=setup_stage5();
iota=cfg.iota;

fprintf('==============================================================\n');
fprintf('Stage 5B: lambda-robustness check on N=8 3-D Morley pullout\n');
fprintf('Pure MATLAB mesh N = 8\n');
fprintf('JS scheme unchanged: P0 stress in J, P0 divergence in s.\n');
fprintf('Corrected nonhomogeneous boundary RHS J(v,g) is retained.\n');
fprintf('mu=%g, iota=%g, lambda = ',cfg.mu,iota);fprintf('%g ',cfg.lambdaList);fprintf('\n');
fprintf('Exact 1-D anti-plane reference is independent of lambda because div u = 0.\n');
fprintf('N=8 uses scaled GMRES with the validated 3-component block preconditioner.\n');
fprintf('The lambda=1 solution is reused only as the initial guess for lambda=1e5.\n');
fprintf('==============================================================\n\n');

pref=stage5_reference_config(cfg);
ref=solve_aifantis_reference_1d(pref,iota);
fprintf('1-D reference max BC residual = %.3e\n',ref.maxBoundaryResidual);
fprintf('inner/outer d_n u3 = %.9g / %.9g\n\n',ref.trace.inner.dnu3,ref.trace.outer.dnu3);

N=8;
m=build_stage4_mesh(N,cfg);
fprintf('Mesh: tet=%d, DOF=%d, free=%d, hmax=%.9g, first dr=%.9g\n', ...
    size(m.elem,1),m.vectorNdof,m.nFree,m.hmax,m.firstRadialCell);
fprintf('Geometry: inner/outer chord mismatch %.3e / %.3e, volume diff %.3e\n', ...
    m.boundary.maxRadiusMismatchInner,m.boundary.maxRadiusMismatchOuter,m.boundary.relativeVolumeDifference);

[gFull,bd]=build_pullout_boundary_dofs(m,ref,cfg);
fprintf('Boundary d_n interpolation rel diff: inner %.3e, outer %.3e\n', ...
    bd.innerFaceMeanRelDiff,bd.outerFaceMeanRelDiff);

baselineFile=fullfile(cfg.resultDir,'stage5a_verified_N4_baseline.csv');
if ~exist(baselineFile,'file')
    error('stage5b:baseline','Missing verified Stage-5A baseline file: %s',baselineFile);
end
B4=readtable(baselineFile);

nLam=numel(cfg.lambdaList);
lambdaCol=zeros(nLam,1); hcol=repmat(m.hmax,nLam,1);
ndof=repmat(m.vectorNdof,nLam,1); nfree=repmat(m.nFree,nLam,1);
L2=zeros(nLam,1); strain=zeros(nLam,1); gradStrain=zeros(nLam,1); threebar=zeros(nLam,1);
H1=zeros(nLam,1); H2=zeros(nLam,1); divUh=zeros(nLam,1); transverse=zeros(nLam,1);
relL2=zeros(nLam,1); relStrain=zeros(nLam,1); relGradStrain=zeros(nLam,1); relThreebar=zeros(nLam,1);
ratioL2=nan(nLam,1); ratioStrain=nan(nLam,1); ratioGradStrain=nan(nLam,1); ratioThreebar=nan(nLam,1);
profileMax=zeros(nLam,1); profileDiff=zeros(nLam,1); profileDiffRel=zeros(nLam,1);
rateL2=nan(nLam,1); rateStrain=nan(nLam,1); rateGradStrain=nan(nLam,1); rateThreebar=nan(nLam,1);
relres=zeros(nLam,1); assemblySec=zeros(nLam,1); solveSec=zeros(nLam,1); errorSec=zeros(nLam,1);
boundaryRhsNorm=zeros(nLam,1); solverName=strings(nLam,1);

baseX=[]; basePath=[]; baseErr=[];
baseProfileScale=max(abs(ref.w));
if baseProfileScale==0, baseProfileScale=cfg.ub; end

for iL=1:nLam
    lambda=cfg.lambdaList(iL);
    lambdaCol(iL)=lambda;
    fprintf('\n================ N=8, lambda=%g ================\n',lambda);

    fprintf('Assembling JS system ...\n');
    ta=tic;
    [A,b,vs]=assemble_volume_js3d_inhom(m,lambda,iota,gFull,cfg);
    [A,b,fs]=assemble_faces_js3d_p0_inhom(m,lambda,gFull,ref,cfg,A,b);
    assemblySec(iL)=toc(ta);
    boundaryRhsNorm(iL)=fs.boundaryDataRhsNorm;
    fprintf('Assembly %.2f s: size=%d, nnz=%d, ||b||=%.6e\n', ...
        assemblySec(iL),size(A,1),nnz(A),norm(b));
    fprintf('Checks: volume sym %.3e, J skew %.3e, S sym %.3e\n', ...
        vs.symmetryDefect,fs.maxJSkewDefect,fs.maxSSymmetryDefect);
    fprintf('Nonhomogeneous boundary J(v,g) RHS norm = %.6e\n',fs.boundaryDataRhsNorm);

    if iL==1
        x0=zeros(m.nFree,1);
    else
        x0=baseX;
        fprintf('Using converged lambda=%g free solution as GMRES initial guess.\n',cfg.lambdaList(1));
    end

    fprintf('Solving ...\n');
    sol=solve_js_system(A,b,cfg,m.nFreeScalar,x0);
    solveSec(iL)=sol.seconds; relres(iL)=sol.relres; solverName(iL)=string(sol.solver);
    fprintf('Solver: %s, time %.2f s, true relative residual %.3e, accepted=%d\n', ...
        sol.solver,sol.seconds,sol.relres,sol.accepted);
    if ~sol.accepted
        error('stage5b:solve','N=8 lambda=%g failed residual acceptance.',lambda);
    end

    uFull=gFull; uFull(double(m.freeDof))=sol.x;
    fixedMismatch=max(abs(uFull(double(m.fixedDof))-gFull(double(m.fixedDof))));
    fprintf('Prescribed-DOF mismatch = %.3e\n',fixedMismatch);
    clear A b x0

    fprintf('Computing approximation errors ...\n');
    te=tic; err=compute_pullout_errors(m,uFull,ref,iota,cfg); errorSec(iL)=toc(te);
    fprintf('Errors: L2=%.6e, strain=%.6e, gradStrain=%.6e, threebar=%.6e\n', ...
        err.L2,err.strain,err.gradStrain,err.threebar);
    fprintf('        H1=%.6e, H2=%.6e, ||div uh||=%.6e\n',err.H1,err.H2,err.divUh);
    fprintf('Relative: L2=%.6e, strain=%.6e, gradStrain=%.6e, threebar=%.6e\n', ...
        err.relL2,err.relStrain,err.relGradStrain,err.relThreebar);
    fprintf('Geometry crescent quadrature-volume fraction = %.3e\n', ...
        err.geometryOutsideReferenceVolumeFraction);

    path=compute_radial_profile(m,uFull,ref,cfg);
    fprintf('Interior radial profile max |u3_h-w| = %.6e, missing points=%d\n', ...
        path.maxAbsError,path.nMissing);

    if iL==1
        baseX=sol.x;
        basePath=path;
        baseErr=err;
        ratioL2(iL)=1; ratioStrain(iL)=1; ratioGradStrain(iL)=1; ratioThreebar(iL)=1;
        profileDiff(iL)=0; profileDiffRel(iL)=0;
    else
        ratioL2(iL)=err.L2/baseErr.L2;
        ratioStrain(iL)=err.strain/baseErr.strain;
        ratioGradStrain(iL)=err.gradStrain/baseErr.gradStrain;
        ratioThreebar(iL)=err.threebar/baseErr.threebar;
        profileDiff(iL)=max(abs(path.uh3-basePath.uh3));
        profileDiffRel(iL)=profileDiff(iL)/baseProfileScale;
        fprintf('Lambda robustness ratios vs lambda=%g on SAME N=8 mesh:\n',cfg.lambdaList(1));
        fprintf('  error ratios: L2=%.4f, strain=%.4f, gradStrain=%.4f, threebar=%.4f\n', ...
            ratioL2(iL),ratioStrain(iL),ratioGradStrain(iL),ratioThreebar(iL));
        fprintf('  max |u3_h(lambda)-u3_h(lambda=1)| = %.6e (relative %.3e)\n', ...
            profileDiff(iL),profileDiffRel(iL));
    end

    % N=4 -> N=8 observed rate using the verified Stage-5A N=4 row.
    ib=find(abs(B4.lambda-lambda)<0.5,1);
    if isempty(ib)
        error('stage5b:baselineLambda','No N=4 baseline row for lambda=%g.',lambda);
    end
    dh=log(B4.hmax(ib)/m.hmax);
    rateL2(iL)=log(B4.L2(ib)/err.L2)/dh;
    rateStrain(iL)=log(B4.strain(ib)/err.strain)/dh;
    rateGradStrain(iL)=log(B4.gradStrain(ib)/err.gradStrain)/dh;
    rateThreebar(iL)=log(B4.threebar(ib)/err.threebar)/dh;
    fprintf('Observed N=4 -> N=8 rate for lambda=%g: L2=%.3f, strain=%.3f, gradStrain=%.3f, threebar=%.3f\n', ...
        lambda,rateL2(iL),rateStrain(iL),rateGradStrain(iL),rateThreebar(iL));

    L2(iL)=err.L2; strain(iL)=err.strain; gradStrain(iL)=err.gradStrain; threebar(iL)=err.threebar;
    H1(iL)=err.H1; H2(iL)=err.H2; divUh(iL)=err.divUh; transverse(iL)=err.transverseL2;
    relL2(iL)=err.relL2; relStrain(iL)=err.relStrain; relGradStrain(iL)=err.relGradStrain; relThreebar(iL)=err.relThreebar;
    profileMax(iL)=path.maxAbsError;

    compact=struct('N',N,'lambda',lambda,'iota',iota,'hmax',m.hmax, ...
        'Ndof',m.vectorNdof,'Nfree',m.nFree,'uFull',uFull,'mesh',m, ...
        'error',err,'path',path,'boundaryDiagnostics',bd,'solver',sol, ...
        'boundaryRhsNorm',fs.boundaryDataRhsNorm, ...
        'assemblySeconds',assemblySec(iL),'errorSeconds',errorSec(iL));
    save(fullfile(cfg.outputDir,sprintf('stage5b_N8_lambda_%g_iota_%g.mat',lambda,iota)), ...
        'compact','-v7.3');
    clear uFull compact err path fs vs sol
end

summary=table(repmat(N,nLam,1),lambdaCol,hcol,ndof,nfree,L2,strain,gradStrain,threebar,H1,H2, ...
    relL2,relStrain,relGradStrain,relThreebar,rateL2,rateStrain,rateGradStrain,rateThreebar, ...
    ratioL2,ratioStrain,ratioGradStrain,ratioThreebar,profileMax,profileDiff,profileDiffRel, ...
    divUh,transverse,repmat(bd.innerFaceMeanRelDiff,nLam,1),repmat(bd.outerFaceMeanRelDiff,nLam,1), ...
    boundaryRhsNorm,relres,assemblySec,solveSec,errorSec,solverName, ...
    'VariableNames',{'N','lambda','hmax','Ndof','Nfree','L2','strain','gradStrain','threebar','H1','H2', ...
    'relL2','relStrain','relGradStrain','relThreebar','rateN4toN8_L2','rateN4toN8_strain', ...
    'rateN4toN8_gradStrain','rateN4toN8_threebar','ratioL2vsLambda1','ratioStrainVsLambda1', ...
    'ratioGradStrainVsLambda1','ratioThreebarVsLambda1','profileMaxError','profileDiffVsLambda1', ...
    'profileDiffRelVsLambda1','divUh','transverseL2','innerNormalRelDiff','outerNormalRelDiff', ...
    'boundaryRhsNorm','relativeResidual','assemblySeconds','solveSeconds','errorSeconds','solver'});

fprintf('\n==============================================================\n');
fprintf('Stage 5B N=8 lambda-robustness summary\n');
disp(summary);

if nLam>=2
    fprintf('N=8 lambda=%g / lambda=%g error ratios: L2=%.4f, strain=%.4f, gradStrain=%.4f, threebar=%.4f\n', ...
        cfg.lambdaList(2),cfg.lambdaList(1),ratioL2(2),ratioStrain(2),ratioGradStrain(2),ratioThreebar(2));
    fprintf('N=8 radial-profile difference relative to max reference displacement: %.3e\n',profileDiffRel(2));
end

csvFile=fullfile(cfg.resultDir,'stage5b_lambda_N8.csv');
matFile=fullfile(cfg.resultDir,'stage5b_lambda_N8.mat');
writetable(summary,csvFile);save(matFile,'summary','cfg','ref','B4');
fprintf('Saved summary CSV: %s\n',csvFile);
fprintf('Saved summary MAT: %s\n',matFile);

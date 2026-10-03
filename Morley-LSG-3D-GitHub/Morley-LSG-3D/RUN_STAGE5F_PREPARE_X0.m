clear; clc;
cfg=setup_stage5();
iota=cfg.iota;
lambdaLo=1;

fprintf('==============================================================\n');
fprintf('Stage 5F-A: prepare trusted N=8 lambda=1 initial guess\n');
fprintf('This stage does NOT alter the finite element formulation.\n');
fprintf('It saves only the N=8 lambda=1 free solution and baseline errors.\n');
fprintf('After this script finishes, CLOSE MATLAB and start a fresh MATLAB\n');
fprintf('session before running RUN_STAGE5F_N8_HIGHLAMBDA.\n');
fprintf('==============================================================\n\n');

pref=stage5_reference_config(cfg);
ref=solve_aifantis_reference_1d(pref,iota);
N=8;
m8=build_stage4_mesh(N,cfg);
[g8,~]=build_pullout_boundary_dofs(m8,ref,cfg);

fprintf('Mesh N=8: tet=%d, DOF=%d, free=%d, hmax=%.9g\n', ...
    size(m8.elem,1),m8.vectorNdof,m8.nFree,m8.hmax);
fprintf('\n--- N=8 lambda=1 baseline solve ---\n');
[A,b]=assemble_volume_js3d_inhom(m8,lambdaLo,iota,g8,cfg);
[A,b]=assemble_faces_js3d_p0_inhom(m8,lambdaLo,g8,ref,cfg,A,b);
cfgLo=cfg;
cfgLo.solver='gmres-block';
cfgLo.solverTol=1e-8;
cfgLo.gmresTol=1e-9;
cfgLo.gmresRestart=40;
cfgLo.gmresMaxit=200;
cfgLo.refinementSteps=0;
sol=solve_js_system(A,b,cfgLo,m8.nFreeScalar,zeros(m8.nFree,1));
if ~sol.accepted
    error('stage5f:x0','N=8 lambda=1 baseline solve failed.');
end
u=g8;
u(double(m8.freeDof))=sol.x;
err8lo=compute_pullout_errors(m8,u,ref,iota,cfg);
path8lo=compute_radial_profile(m8,u,ref,cfg);
x0=sol.x; %#ok<NASGU>
hmax8=m8.hmax; %#ok<NASGU>
nFree8=m8.nFree; %#ok<NASGU>

fprintf('N=8 lambda=1 errors: L2=%.6e, strain=%.6e, gradStrain=%.6e, threebar=%.6e\n', ...
    err8lo.L2,err8lo.strain,err8lo.gradStrain,err8lo.threebar);
fprintf('true relative residual = %.3e\n',sol.relres);

saveFile=fullfile(cfg.resultDir,'stage5f_N8_lambda1_x0.mat');
save(saveFile,'x0','err8lo','path8lo','hmax8','nFree8','-v7.3');
fprintf('\nSaved initial guess: %s\n',saveFile);
fprintf('Stage 5F-A complete. CLOSE MATLAB now, then start a fresh MATLAB\n');
fprintf('session and run RUN_STAGE5F_N8_HIGHLAMBDA.\n');

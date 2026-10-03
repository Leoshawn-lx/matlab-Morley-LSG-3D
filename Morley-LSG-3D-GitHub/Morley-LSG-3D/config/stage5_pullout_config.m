function cfg = stage5_pullout_config()
%STAGE5_PULLOUT_CONFIG Solver-accuracy diagnostic for lambda=1e5.
%
% This stage does NOT change the finite element formulation.  It only tightens
% the N=8 Krylov solve and adds residual-correction iterative refinement.

cfg.root = fileparts(fileparts(mfilename('fullpath')));
cfg.codeVersion = 'morley-pullout-stage5f-fresh-N8-exactH';

% Geometry/loading.
cfg.RI = 0.01;
cfg.RE = 1.00;
cfg.H  = 0.05;
cfg.ub = 0.005;

% Parameters.
cfg.mu = 1;
cfg.iota = 0.1;
cfg.lambdaList = [1,1e5];
% Compatibility list used by the Stage-1 reference adapter.  Stage 5C
% actually solves only cfg.iota, but the reference helper expects iotaList.
cfg.iotaList = cfg.iota;

% Mesh hierarchy used only by the diagnostic runner.
cfg.NList  = [4,8];
% Structured annular mesh parameters corresponding to N=4 and N=8.
cfg.nrList = [4,8];
cfg.ntList = [24,48];
cfg.nzList = [2,4];
cfg.maxFreeDof = 500000;
cfg.warnFreeDof = 150000;

% Validated Stage-1 reference settings.
cfg.reference.nInit = 240;
cfg.reference.nPlot = 6001;
cfg.reference.RelTol = 1e-9;
cfg.reference.AbsTol = 1e-11;
cfg.reference.NMax = 60000;

% Assembly/boundary quadrature.
cfg.edgeQuadOrder = 8;
cfg.faceQuadOrder = 6;
cfg.chunkElements = 128;
cfg.chunkFaces = 256;

% Error/postprocessing.
cfg.errorOrder = 6;
cfg.profilePoints = 301;

% Stronger N=8 solver diagnostic.
cfg.solver = 'gmres-block';
cfg.solverTol = 1e-12;
cfg.gmresTol = 1e-10;
cfg.gmresRestart = 40;
cfg.gmresMaxit = 300;
cfg.blockIcholDropTolList = [1e-2,5e-3];
cfg.blockIcholDiagCompList = [0,1e-3,1e-2,1e-1,1,10];

% Residual-correction iterative refinement.  Each correction solves
%   A*delta = b-A*x
% with the SAME preconditioner; A itself is never altered.
cfg.refinementSteps = 2;
cfg.refinementTol = 1e-10;
cfg.refinementMaxit = 300;

cfg.outputDir = fullfile(cfg.root,'output');
cfg.resultDir = fullfile(cfg.root,'results');
end

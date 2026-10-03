function cfg=setup_stage5()
root=fileparts(mfilename('fullpath'));
addpath(root);
addpath(fullfile(root,'config'));
addpath(fullfile(root,'mesh'));
addpath(fullfile(root,'fem'));
addpath(fullfile(root,'reference'));
addpath(fullfile(root,'assembly'));
addpath(fullfile(root,'solver'));
addpath(fullfile(root,'post'));
addpath(fullfile(root,'utils'));
cfg=stage5_pullout_config();
if ~exist(cfg.outputDir,'dir'),mkdir(cfg.outputDir);end
if ~exist(cfg.resultDir,'dir'),mkdir(cfg.resultDir);end
end

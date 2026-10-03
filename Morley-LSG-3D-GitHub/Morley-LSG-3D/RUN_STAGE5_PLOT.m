clear; clc;
cfg=setup_stage5();
matFile=fullfile(cfg.resultDir,'stage5b_lambda_N8.mat');
if ~exist(matFile,'file'), error('Run RUN_STAGE5_LAMBDA_N8 first.'); end
S=load(matFile,'summary'); T=S.summary;

figure('Color','w','Name','N=8 lambda robustness errors');
semilogy(T.lambda,T.L2,'o-','LineWidth',1.4); hold on;
semilogy(T.lambda,T.strain,'s-','LineWidth',1.4);
semilogy(T.lambda,T.gradStrain,'^-','LineWidth',1.4);
semilogy(T.lambda,T.threebar,'d-','LineWidth',1.4);
xlabel('\lambda'); ylabel('absolute error'); grid on; box on;
legend('L2','strain','grad-strain','three-bar','Location','best');
set(gca,'XScale','log');

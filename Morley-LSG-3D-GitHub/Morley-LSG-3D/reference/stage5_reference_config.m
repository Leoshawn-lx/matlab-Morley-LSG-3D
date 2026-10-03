function p = stage5_reference_config(cfg)
% Adapter to the validated Stage-1 reference routine.
p.RI=cfg.RI; p.RE=cfg.RE; p.H=cfg.H; p.ub=cfg.ub;
p.mu=cfg.mu; p.lambdaList=cfg.lambdaList;
if isfield(cfg,'iotaList')
    p.iotaList=cfg.iotaList;
else
    p.iotaList=cfg.iota;
end
p.reference=cfg.reference;
end

function s=cfg_signature(cfg)
% Minimal deterministic signature consumed by the Stage-2 topology helper.
s=sprintf('RI=%.17g|RE=%.17g|H=%.17g|ub=%.17g|mu=%.17g|ver=%s', ...
    cfg.RI,cfg.RE,cfg.H,cfg.ub,cfg.mu,cfg.codeVersion);
end

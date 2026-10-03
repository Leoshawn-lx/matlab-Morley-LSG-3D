function m = build_stage4_mesh(N,cfg)
%BUILD_STAGE4_MESH Generate and preprocess one pure-MATLAB pullout mesh.
idx = find(cfg.NList==N,1);
if isempty(idx), error('stage4:meshLevel','Unsupported N=%g.',N); end
nr=cfg.nrList(idx); nt=cfg.ntList(idx); nz=cfg.nzList(idx);
raw=build_annular_structured_tet_mesh(nr,nt,nz,cfg);
m=build_morley3d_topology(raw,cfg);
m.N=N; m.nr=nr; m.nt=nt; m.nz=nz;
m.radii=raw.radii; m.theta=raw.theta; m.zlev=raw.zlev;
m.firstRadialCell=raw.firstRadialCell;
m.boundary=pullout_mesh_diagnostics(m,cfg);
end

function path = compute_radial_profile(m,uFull,ref,cfg)
%COMPUTE_RADIAL_PROFILE Interior radial displacement profile.
% The path is close to theta=0 and z=H/2, but nudged away from mesh faces
% so that the nonconforming polynomial trace is unambiguous.

npt=cfg.profilePoints;
r=linspace(cfg.RI,cfg.RE,npt)';
theta=1e-10;
z=(0.5+1e-9)*cfg.H;
pts=[r*cos(theta),r*sin(theta),z*ones(npt,1)];
epsr=1e-9*(cfg.RE-cfg.RI);
pts(1,1:2)=(cfg.RI+epsr)*[cos(theta),sin(theta)];
pts(end,1:2)=(cfg.RE-epsr)*[cos(theta),sin(theta)];
TR=triangulation(double(m.elem),m.node);
tid=pointLocation(TR,pts);
U=nan(npt,3);
for q=1:npt
    if isnan(tid(q)),continue;end
    t=tid(q);
    [~,fd]=element_dofs(m,t);
    coeff=reshape(uFull(fd(:)),10,3);
    phi=morley_basis3d(pts(q,:),element_geometry(m,t));
    U(q,:)=phi*coeff;
end
[w,~,~,~]=eval_aifantis_reference_1d(ref,r);
U(1,:)=[0,0,cfg.ub];U(end,:)=[0,0,0];
path.r=r;path.uh=U;path.uh3=U(:,3);path.uref3=w;
d=abs(path.uh3-path.uref3); d=d(isfinite(d));
if isempty(d),path.maxAbsError=NaN;else,path.maxAbsError=max(d);end
path.nMissing=nnz(isnan(path.uh3));
end

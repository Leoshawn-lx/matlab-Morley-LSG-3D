function ex = pullout_reference_fields(x,ref)
%PULLOUT_REFERENCE_FIELDS Exact radial displacement, gradient and Hessian.
%
% Array conventions:
%   ex.u(q,c)          = u_c(x_q)
%   ex.grad(q,c,j)     = d_j u_c(x_q)
%   ex.hess(q,c,j,k)   = d_j d_k u_c(x_q)
%
% The straight-sided inner polygon contains a very small crescent with
% r<RI. As in Stage 3, radial evaluation is clamped to [RI,RE]. The size of
% that geometry-mismatch region is reported separately by the error routine.

rActual=hypot(x(:,1),x(:,2));
rEval=min(max(rActual,ref.r(1)),ref.r(end));
[w,wp,wpp,~]=eval_aifantis_reference_1d(ref,rEval);

% Radial direction uses the physical point direction. rActual is positive
% throughout the annular mesh.
erx=x(:,1)./rActual;
ery=x(:,2)./rActual;

wx=wp.*erx;
wy=wp.*ery;
wxx=wpp.*erx.^2 + (wp./rEval).*(1-erx.^2);
wyy=wpp.*ery.^2 + (wp./rEval).*(1-ery.^2);
wxy=(wpp-wp./rEval).*erx.*ery;

n=size(x,1);
ex.u=zeros(n,3);
ex.u(:,3)=w;
ex.grad=zeros(n,3,3);
ex.grad(:,3,1)=wx;
ex.grad(:,3,2)=wy;
ex.hess=zeros(n,3,3,3);
ex.hess(:,3,1,1)=wxx;
ex.hess(:,3,1,2)=wxy;
ex.hess(:,3,2,1)=wxy;
ex.hess(:,3,2,2)=wyy;
ex.rActual=rActual;
ex.rEval=rEval;
end

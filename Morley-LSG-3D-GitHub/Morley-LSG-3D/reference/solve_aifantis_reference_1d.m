function ref = solve_aifantis_reference_1d(p,iota)
%SOLVE_AIFANTIS_REFERENCE_1D Solve the radial anti-plane Aifantis reference.
%
%   u_ref(x,y,z) = (0,0,w(r))^T,  r = sqrt(x^2+y^2).
%
% The scalar equation is
%
%   iota^2 Delta_r^2 w - Delta_r w = 0,
%   Delta_r w = w'' + w'/r,
%
% with
%
%   w(RI)  = ub,      w(RE)  = 0,
%   w''(RI)= 0,       w''(RE)= 0.
%
% State vector y = [w,w',w'',w''']^T.  The fourth-order ODE is
%
%   w'''' = (w''+w'/r)/iota^2 - 2 w'''/r
%           + w''/r^2 - w'/r^3.
%
% The returned structure contains the bvp4c solution, dense sampled data,
% interpolants, and the nonhomogeneous essential boundary traces that will
% be imposed in the later 3-D Morley calculation.

if ~(isscalar(iota) && isfinite(iota) && iota > 0)
    error('iota must be a positive finite scalar.');
end

RI = p.RI;
RE = p.RE;
ub = p.ub;

odefun = @(r,y) [ ...
    y(2); ...
    y(3); ...
    y(4); ...
    (y(3) + y(2)./r)./iota^2 - 2*y(4)./r + y(3)./r.^2 - y(2)./r.^3];

bcfun = @(ya,yb) [ya(1)-ub; yb(1); ya(3); yb(3)];

% Local-elasticity solution is a good interior initial guess.  It does not
% satisfy w''=0 at the endpoints, which is exactly the boundary-layer
% correction that bvp4c resolves.
L = log(RE/RI);
guess = @(r) [ ...
    ub*log(RE./r)/L; ...
   -ub./(r*L); ...
    ub./(r.^2*L); ...
   -2*ub./(r.^3*L)];

% Cosine clustering resolves both radial endpoints without hard-coding a
% length-scale-dependent mesh.
s = linspace(0,1,p.reference.nInit);
rInit = RI + (RE-RI)*0.5*(1-cos(pi*s));
solinit = bvpinit(rInit,guess);
opts = bvpset('RelTol',p.reference.RelTol, ...
              'AbsTol',p.reference.AbsTol, ...
              'NMax',p.reference.NMax, ...
              'Stats','off');
sol = bvp4c(odefun,bcfun,solinit,opts);

sDense = linspace(0,1,p.reference.nPlot);
r = RI + (RE-RI)*0.5*(1-cos(pi*sDense));
y = deval(sol,r);

ref.iota = iota;
ref.sol = sol;
ref.r = r(:);
ref.w = y(1,:).';
ref.wp = y(2,:).';
ref.wpp = y(3,:).';
ref.wppp = y(4,:).';

ref.wInterp    = griddedInterpolant(ref.r,ref.w,   'pchip','nearest');
ref.wpInterp   = griddedInterpolant(ref.r,ref.wp,  'pchip','nearest');
ref.wppInterp  = griddedInterpolant(ref.r,ref.wpp, 'pchip','nearest');
ref.wpppInterp = griddedInterpolant(ref.r,ref.wppp,'pchip','nearest');

% Boundary residuals of the defining 1-D BVP.
ref.boundaryResidual = [ ...
    ref.w(1)-ub; ...
    ref.w(end); ...
    ref.wpp(1); ...
    ref.wpp(end)];

% Exact essential traces to be used by the 3-D Morley lifting.
% Inner-cylinder outward normal of the annulus is -e_r.
ref.trace.inner.u3  = ref.w(1);
ref.trace.inner.dnu3 = -ref.wp(1);
% Outer-cylinder outward normal is +e_r.
ref.trace.outer.u3  = ref.w(end);
ref.trace.outer.dnu3 = ref.wp(end);
% On z=0,H the reference is independent of z.
ref.trace.plane.dnu3 = 0;

ref.maxBoundaryResidual = max(abs(ref.boundaryResidual));
if ref.maxBoundaryResidual > 1e-7
    warning('iota=%g: 1-D boundary residual is %.3e.', ...
        iota,ref.maxBoundaryResidual);
end

% Simple qualitative checks.  The radial displacement should remain within
% the imposed displacement range up to a small numerical tolerance.
rangeTol = 1e-8*max(1,abs(ub));
ref.rangeViolation = max([max(ref.w-ub), max(-ref.w), 0]);
if ref.rangeViolation > rangeTol
    warning('iota=%g: w(r) leaves [0,ub] by %.3e.', ...
        iota,ref.rangeViolation);
end
end

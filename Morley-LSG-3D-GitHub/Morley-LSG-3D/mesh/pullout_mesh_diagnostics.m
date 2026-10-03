function d = pullout_mesh_diagnostics(m,cfg)
% Geometry/topology diagnostics only. No manufactured solution is used.
% Boundary type convention inherited from build_morley3d_topology:
%   1=Inner, 2=Outer, 3=Bottom, 4=Top.

[b,w] = local_triangle_rule();
d = struct();
d.maxRadiusMismatchInner = 0;
d.maxRadiusMismatchOuter = 0;
d.maxZMistmatchBottom = 0;
d.maxZMistmatchTop = 0;
d.minOutwardDotInner = inf;
d.minOutwardDotOuter = inf;
d.minOutwardDotBottom = inf;
d.minOutwardDotTop = inf;
d.nFaceByGroup = zeros(1,4);

for kk = 1:numel(m.boundaryFaces)
    f = double(m.boundaryFaces(kk));
    typ = double(m.boundaryType(f));
    d.nFaceByGroup(typ) = d.nFaceByGroup(typ) + 1;
    p = m.node(double(m.face(f,:)),:);
    x = b*p;
    n = m.faceNormal(f,:);
    switch typ
        case 1 % inner cylinder; outward domain normal points toward the axis
            rr = sqrt(sum(x(:,1:2).^2,2));
            d.maxRadiusMismatchInner = max(d.maxRadiusMismatchInner,max(abs(rr-cfg.RI)));
            er = [mean(x(:,1)),mean(x(:,2)),0]; er = er/max(norm(er),realmin);
            d.minOutwardDotInner = min(d.minOutwardDotInner,dot(n,-er));
        case 2 % outer cylinder
            rr = sqrt(sum(x(:,1:2).^2,2));
            d.maxRadiusMismatchOuter = max(d.maxRadiusMismatchOuter,max(abs(rr-cfg.RE)));
            er = [mean(x(:,1)),mean(x(:,2)),0]; er = er/max(norm(er),realmin);
            d.minOutwardDotOuter = min(d.minOutwardDotOuter,dot(n,er));
        case 3 % bottom
            d.maxZMistmatchBottom = max(d.maxZMistmatchBottom,max(abs(x(:,3))));
            d.minOutwardDotBottom = min(d.minOutwardDotBottom,dot(n,[0 0 -1]));
        case 4 % top
            d.maxZMistmatchTop = max(d.maxZMistmatchTop,max(abs(x(:,3)-cfg.H)));
            d.minOutwardDotTop = min(d.minOutwardDotTop,dot(n,[0 0 1]));
    end
end

exactVol = pi*(cfg.RE^2-cfg.RI^2)*cfg.H;
d.relativeVolumeDifference = abs(m.totalVolume-exactVol)/exactVol;
d.nBoundaryFaces = numel(m.boundaryFaces);
d.nBoundaryEdges = numel(m.boundaryEdges);
d.scalarBoundaryDof = numel(m.boundaryEdges)+numel(m.boundaryFaces);
d.vectorBoundaryDof = 3*d.scalarBoundaryDof;
d.orientationOK = min([d.minOutwardDotInner,d.minOutwardDotOuter, ...
    d.minOutwardDotBottom,d.minOutwardDotTop]) > 0;

% The quadrature weights are deliberately referenced so MATLAB does not
% warn about an unused output if this routine is modified later.
d.quadWeightSum = sum(w);
end

function [b,w] = local_triangle_rule()
% Seven-point symmetric rule (barycentric coordinates), sufficient for
% geometry sampling; weights sum to one.
a1 = 1/3;
a2 = 0.059715871789770; b2 = 0.470142064105115;
a3 = 0.797426985353087; b3 = 0.101286507323456;
b = [a1 a1 a1; ...
     a2 b2 b2; b2 a2 b2; b2 b2 a2; ...
     a3 b3 b3; b3 a3 b3; b3 b3 a3];
w = [0.225; repmat(0.132394152788506,3,1); repmat(0.125939180544827,3,1)];
end

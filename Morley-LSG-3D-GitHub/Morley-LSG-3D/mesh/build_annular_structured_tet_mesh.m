function raw = build_annular_structured_tet_mesh(nr,nt,nz,cfg)
%BUILD_ANNULAR_STRUCTURED_TET_MESH Pure-MATLAB annular-cylinder tetra mesh.
%
% The construction follows the previously validated pullout hierarchy:
%  1) logarithmic radial rings;
%  2) periodic angular triangulation with alternating diagonals;
%  3) extrusion through z;
%  4) each triangular prism is split into three tetrahedra.
%
% Output deliberately mimics the small subset of the old Gmsh reader
% structure consumed by build_morley3d_topology:
%   raw.node, raw.tet, raw.tri, raw.triTag, raw.tetTag,
%   raw.physicalNames.
%
% Physical tags:
%   101 Inner, 102 Outer, 103 Bottom, 104 Top, 201 Domain.

validateattributes(nr,{'numeric'},{'scalar','integer','positive'});
validateattributes(nt,{'numeric'},{'scalar','integer','>=',4});
validateattributes(nz,{'numeric'},{'scalar','integer','positive'});
if cfg.RI <= 0 || cfg.RE <= cfg.RI || cfg.H <= 0
    error('pullout:geometry','Require 0 < RI < RE and H > 0.');
end

% Logarithmic radial spacing: exactly the spacing used by the earlier
% pullout code and particularly useful near RI=0.01.
radii = cfg.RI*(cfg.RE/cfg.RI).^((0:nr)/nr);
theta = 2*pi*(0:nt-1)/nt;

% One 2-D annular layer.
n2 = (nr+1)*nt;
P = zeros(n2,2);
idx = @(ir,it) ir*nt + it;  % ir=0:nr, it=1:nt
for ir = 0:nr
    for it = 1:nt
        id = idx(ir,it);
        P(id,:) = radii(ir+1)*[cos(theta(it)),sin(theta(it))];
    end
end

% Triangulate each annular quadrilateral. Alternating diagonals avoid a
% strong one-direction bias while remaining periodic at theta=0/2*pi.
tri2 = zeros(2*nr*nt,3);
q = 0;
for ir = 1:nr
    for it = 1:nt
        jp = mod(it,nt)+1;
        a = idx(ir-1,it); b = idx(ir-1,jp);
        c = idx(ir,it);   d = idx(ir,jp);
        if mod(ir+it,2)==0
            pair = [a,c,d; a,d,b];
        else
            pair = [a,c,b; b,c,d];
        end
        for s = 1:2
            one = pair(s,:);
            X = P(one,:);
            if det([X(2,:)-X(1,:);X(3,:)-X(1,:)]) < 0
                one([2,3]) = one([3,2]);
            end
            q = q+1;
            tri2(q,:) = one;
        end
    end
end

% Extrude the layer through the thickness.
zlev = linspace(0,cfg.H,nz+1);
node = zeros((nz+1)*n2,3);
for iz = 0:nz
    ids = iz*n2+(1:n2);
    node(ids,:) = [P,zlev(iz+1)*ones(n2,1)];
end

% Split every triangular prism into three tetrahedra. The sorted base
% triangle convention is inherited from the earlier validated pullout code.
tet = zeros(nz*size(tri2,1)*3,4);
e = 0;
for iz = 1:nz
    off = (iz-1)*n2;
    for kt = 1:size(tri2,1)
        v = sort(tri2(kt,:));
        bot = v+off;
        top = bot+n2;
        prism = [bot(1),bot(2),bot(3),top(3); ...
                 bot(1),bot(2),top(2),top(3); ...
                 bot(1),top(1),top(2),top(3)];
        for s = 1:3
            one = prism(s,:);
            X = node(one,:);
            J = [X(2,:)-X(1,:);X(3,:)-X(1,:);X(4,:)-X(1,:)];
            if det(J) < 0
                one([2,3]) = one([3,2]);
            end
            e = e+1;
            tet(e,:) = one;
        end
    end
end

% Extract the actual topological boundary triangles from the tetrahedra.
% This guarantees exact agreement between the boundary list and the volume
% topology, avoiding any independent surface-mesh convention.
lf = [2 3 4;1 4 3;1 2 4;1 3 2];
ntet = size(tet,1);
allF = zeros(4*ntet,3);
for k = 1:4
    allF((k-1)*ntet+(1:ntet),:) = sort(tet(:,lf(k,:)),2);
end
[face,~,ic] = unique(allF,'rows');
counts = accumarray(ic,1,[size(face,1),1]);
bnd = face(counts==1,:);

% Classify each topological boundary triangle. Since all generating nodes
% lie exactly on r=RI/RE or z=0/H, vertex-based classification is exact up
% to roundoff. The curved cylinder itself remains the usual chordal planar
% approximation associated with first-order tetrahedra.
% IMPORTANT: bnd is an nBoundary-by-3 matrix of vertex indices.  In
% MATLAB, node(bnd,1) does not reliably preserve the nBoundary-by-3 face
% layout; the subscript matrix may be linearized.  Rebuild the vertex
% coordinate arrays explicitly so classification is performed face-by-face,
% not vertex-by-vertex.
xb = reshape(node(bnd(:),1),size(bnd));
yb = reshape(node(bnd(:),2),size(bnd));
zb = reshape(node(bnd(:),3),size(bnd));
rb = sqrt(xb.^2 + yb.^2);
tolR = 5e-11*max(1,cfg.RE);
tolZ = 5e-11*max(1,cfg.H);
allInner = all(abs(rb-cfg.RI)<=tolR,2);
allOuter = all(abs(rb-cfg.RE)<=tolR,2);
allBottom = all(abs(zb)<=tolZ,2);
allTop = all(abs(zb-cfg.H)<=tolZ,2);
which = double(allInner) + 2*double(allOuter) + 3*double(allBottom) + 4*double(allTop);
numClass = double(allInner)+double(allOuter)+double(allBottom)+double(allTop);
if any(numClass~=1)
    bad = find(numClass~=1,1);
    error('pullout:boundaryClassification', ...
        ['Boundary face %d = [%d %d %d] was assigned to %d groups. ', ...
         'flags [Inner Outer Bottom Top] = [%d %d %d %d].'], ...
        bad,bnd(bad,1),bnd(bad,2),bnd(bad,3),numClass(bad), ...
        allInner(bad),allOuter(bad),allBottom(bad),allTop(bad));
end
triTag = zeros(size(bnd,1),1);
triTag(allInner)  = 101;
triTag(allOuter)  = 102;
triTag(allBottom) = 103;
triTag(allTop)    = 104;

physicalNames = struct('dim',{},'tag',{},'name',{});
physicalNames(1) = struct('dim',2,'tag',101,'name','Inner');
physicalNames(2) = struct('dim',2,'tag',102,'name','Outer');
physicalNames(3) = struct('dim',2,'tag',103,'name','Bottom');
physicalNames(4) = struct('dim',2,'tag',104,'name','Top');
physicalNames(5) = struct('dim',3,'tag',201,'name','Domain');

raw = struct();
raw.node = node;
raw.tet = tet;
raw.tri = bnd;
raw.triTag = triTag;
raw.tetTag = 201*ones(size(tet,1),1);
raw.physicalNames = physicalNames;
raw.radii = radii(:);
raw.theta = theta(:);
raw.zlev = zlev(:);
raw.nr = nr;
raw.nt = nt;
raw.nz = nz;
raw.firstRadialCell = radii(2)-radii(1);
raw.minRadialCell = min(diff(radii));
raw.maxRadialCell = max(diff(radii));
end

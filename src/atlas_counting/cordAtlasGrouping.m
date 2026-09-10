function grouping = cordAtlasGrouping(av, segmentinfo, nativeNz)
%CORDATLASGROUPING Index spinal cord atlas voxels by region and by segment.
%
%   GROUPING = CORDATLASGROUPING(AV, SEGMENTINFO) builds the two lookups every
%   spinal cord statistic needs: which region a voxel belongs to, and which
%   spinal segment (C1, C2, ... ) it falls in. Regions come from the annotation
%   volume AV, segments from the Start/End planes in SEGMENTINFO
%   (Segments.csv). Cord results are always reported per region *and* per
%   segment, because the same region runs the whole length of the cord and
%   averaging over it would throw away the rostrocaudal axis.
%
%   GROUPING = CORDATLASGROUPING(AV, SEGMENTINFO, NATIVENZ) says that the
%   segment planes in SEGMENTINFO are indices into a volume of NATIVENZ planes,
%   while AV is on a grid of its own. Every plane of AV is mapped back to the
%   native plane it came from and assigned the segment holding it, so statistics
%   can be computed on the registration grid as easily as on the native one.
%   NATIVENZ defaults to size(AV, 3).
%
%   OUTPUT
%     grouping - struct with fields
%        avinds     - Nareas x 1 region ids present in AV, sorted (0 included,
%                     as everything outside the cord).
%        lut        - lookup from region id + 1 to its row in avinds, so
%                     lut(av + 1) gives row indices directly.
%        segofplane - Nz x 1 segment index of every plane of AV, 0 for planes
%                     no segment covers.
%        segnames   - Nsegments x 1 segment names.
%        Nareas, Nsegments, volsize
%
%   The lookups are deliberately returned rather than applied: the consumers
%   walk the volume one segment at a time, which keeps a whole-volume index
%   array off the heap.
%
%   See also CORDAREASTATISTICS, GROUPCORDPOINTSINTOAREAS, LOADSPINALCORDATLAS.

if nargin < 3 || isempty(nativeNz)
    nativeNz = size(av, 3);
end

grouping         = struct();
grouping.volsize = size(av);
grouping.avinds  = unique(av(:));                % one output: no per-voxel index

lut = zeros(double(max(grouping.avinds)) + 1, 1, 'uint32');
lut(double(grouping.avinds) + 1) = 1:numel(grouping.avinds);
grouping.lut    = lut;
grouping.Nareas = numel(grouping.avinds);

%--------------------------------------------------------------------------
% assign every plane of av to a segment, going through the native grid
Nz     = size(av, 3);
zscale = Nz / nativeNz;

% centre of grid plane g sits at native plane (g - 0.5)/zscale + 0.5
nativeplane = round(((1:Nz)' - 0.5) / zscale + 0.5);

segstart = double(segmentinfo.Start(:));
segend   = double(segmentinfo.End(:));
Nsegments = numel(segstart);

segofplane = zeros(Nz, 1);
for ii = 1:Nsegments
    inseg = nativeplane >= segstart(ii) & nativeplane <= segend(ii);
    segofplane(inseg) = ii;
end

grouping.segofplane = segofplane;
grouping.Nsegments  = Nsegments;
if ismember('Segment', segmentinfo.Properties.VariableNames)
    grouping.segnames = string(segmentinfo.Segment);
else
    grouping.segnames = string((1:Nsegments)');
end

end

function [areacounts, areavols, ptareas, ptsegments] = groupCordPointsIntoAreas(...
    cellcoords, av, grouping, atlasres)
%GROUPCORDPOINTSINTOAREAS Count points per spinal cord region and segment.
%
%   [AREACOUNTS, AREAVOLS, PTAREAS, PTSEGMENTS] = GROUPCORDPOINTSINTOAREAS(
%   CELLCOORDS, AV, GROUPING, ATLASRES) counts the points CELLCOORDS - already
%   rounded to atlas voxels, [x y z] - into the regions of the annotation AV and
%   the spinal segments described by GROUPING (see CORDATLASGROUPING).
%
%   This is the cord counterpart of GROUPCELLSINTOLEAFREGIONS. A brain is split
%   into hemispheres there; a cord is split along its length instead, so counts
%   come out as Nareas x Nsegments rather than Nareas x 2.
%
%   AREAVOLS is the volume of every region in every segment, in mm^3, so counts
%   can be turned into densities. PTAREAS and PTSEGMENTS give the region id and
%   the segment index of every input point, for anything that needs to work
%   point by point; a segment index of 0 means the point fell on a plane no
%   segment covers.
%
%   See also CORDATLASGROUPING, CORDAREASTATISTICS, TRANSFORMCORDPOINTSTOATLAS.

Nareas    = grouping.Nareas;
Nsegments = grouping.Nsegments;

%--------------------------------------------------------------------------
% region and segment of every point
idx     = sub2ind(size(av), cellcoords(:, 2), cellcoords(:, 1), cellcoords(:, 3));
ptareas = av(idx);
arearow = double(grouping.lut(double(ptareas) + 1));

ptsegments = grouping.segofplane(cellcoords(:, 3));

%--------------------------------------------------------------------------
ikeep      = arearow > 0 & ptsegments > 0;
areacounts = accumarray([arearow(ikeep) ptsegments(ikeep)], 1, ...
    [Nareas Nsegments], @sum);

%--------------------------------------------------------------------------
% region volumes, walked one segment at a time to keep memory flat
areavols = zeros(Nareas, Nsegments, 'single');
voxvol   = prod(atlasres) * 1e-9;   % um^3 -> mm^3
for iseg = 1:Nsegments
    iplanes = grouping.segofplane == iseg;
    if ~any(iplanes)
        continue
    end
    avblock = av(:, :, iplanes);
    subs    = double(grouping.lut(double(avblock(:)) + 1));
    areavols(:, iseg) = single(accumarray(subs, 1, [Nareas 1], @sum)) * voxvol;
end

end

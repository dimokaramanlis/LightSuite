function [volumecounts, volumeareas, cellatlasids] = groupCellsIntoLeafRegions(cellcoords, annvol, groupinds)
%UNTITLED Summary of this function goes here
%   volumeareas is in mm3

% we assume the annotation volume resolution is 10um per voxel

Ngroups = max(annvol, [], 'all')+1;
%--------------------------------------------------------------------------
ileft          = cellcoords(:, 3) > size(annvol, 3)/2;
linearIndicesL = sub2ind(size(annvol),  cellcoords(ileft, 2),  cellcoords(ileft, 1),  cellcoords(ileft, 3));
linearIndicesR = sub2ind(size(annvol), cellcoords(~ileft, 2), cellcoords(~ileft, 1), cellcoords(~ileft, 3));
cellatlasidsR  = annvol(linearIndicesR);
cellatlasidsL  = annvol(linearIndicesL);
cellatlasids   = {[find(~ileft) uint64(cellatlasidsR)], [find(ileft) uint64(cellatlasidsL)]};
%--------------------------------------------------------------------------
% calculate areas
volumeareas = accumarray(reshape(annvol+1, numel(annvol), 1), 1, [Ngroups 1], @sum);
volumeareas = volumeareas*1e-6; % in mm3
%--------------------------------------------------------------------------
% calculate counts
volumecountsR = accumarray(cellatlasidsR+1, 1, [Ngroups 1], @sum);
volumecountsL = accumarray(cellatlasidsL+1, 1, [Ngroups 1], @sum);
volumecounts  = [volumecountsR volumecountsL];
%--------------------------------------------------------------------------
volumecounts = volumecounts(groupinds+1,:);
volumeareas  = volumeareas(groupinds+1,:);
%--------------------------------------------------------------------------
end
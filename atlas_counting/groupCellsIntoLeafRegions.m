function [volumecounts, volumeareas] = groupCellsIntoLeafRegions(cellcoords, annvol)
%UNTITLED Summary of this function goes here
%   volumeareas is in mm3

% we assume the annotation volume resolution is 10um per voxel

Ngroups = max(annvol, [], 'all');
%--------------------------------------------------------------------------
linearIndices = sub2ind(size(annvol), cellcoords(:, 2), cellcoords(:, 1), cellcoords(:, 3));
cellatlasids  = annvol(linearIndices);
%--------------------------------------------------------------------------
% calculate areas
volumeareas = accumarray(reshape(annvol, numel(annvol), 1), 1, [Ngroups 1], @sum);
volumeareas = volumeareas*1e-6; % in mm3
%--------------------------------------------------------------------------
% calculate counts
volumecounts = accumarray(cellatlasids, 1, [Ngroups 1], @sum);
%--------------------------------------------------------------------------

end
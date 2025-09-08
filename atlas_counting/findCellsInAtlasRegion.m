function [idxarea] = findCellsInAtlasRegion(cellcoords, annvol, groupinds)
%UNTITLED Summary of this function goes here
%   volumeareas is in mm3
%--------------------------------------------------------------------------
cellcoords = round(cellcoords(:,1:3));

irem0 = any(cellcoords(:,1:3)<1, 2);
iremx = cellcoords(:,1) > size(annvol, 2);
iremy = cellcoords(:,2) > size(annvol, 1);
iremz = cellcoords(:,3) > size(annvol, 3);
irem  = irem0 | iremx | iremy | iremz;
ikeep = find(~irem);
%--------------------------------------------------------------------------
cellinds = sub2ind(size(annvol),  cellcoords(ikeep, 2),  cellcoords(ikeep, 1),  cellcoords(ikeep, 3));
iarea = ismembc(annvol(cellinds), uint16(groupinds));
idxarea = ikeep(iarea);
%--------------------------------------------------------------------------
end
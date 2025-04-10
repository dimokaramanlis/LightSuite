function finalpts = sanitizeCellCoords(cellcoords, annvol)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

finalpts = round(cellcoords); % rounding to transfer correspondance to pixels
%--------------------------------------------------------------------------
% remove cells outside volume
irem0 = any(finalpts<1, 2);
iremx = finalpts(:,1) > size(annvol, 2);
iremy = finalpts(:,2) > size(annvol, 1);
iremz = finalpts(:,3) > size(annvol, 3);
irem  = irem0 | iremx | iremy | iremz;
finalpts(irem, :) = [];
%--------------------------------------------------------------------------
% remove cells outside brain regions (root)
indcells = sub2ind(size(annvol), finalpts(:,2), finalpts(:,1), finalpts(:,3));
irem     = annvol(indcells) == 1;
finalpts(irem, :) = [];
%--------------------------------------------------------------------------
end
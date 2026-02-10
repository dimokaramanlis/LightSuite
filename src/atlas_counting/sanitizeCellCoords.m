function finalpts = sanitizeCellCoords(cellcoords, annvol, varargin)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%--------------------------------------------------------------------------
if nargin < 3
    igood = true(size(cellcoords, 1), 1);
else
    igood = varargin{1};
end
%--------------------------------------------------------------------------
finalpts = [round(cellcoords(:, 1:3)) cellcoords(:, 4:end)]; % rounding to transfer correspondance to pixels
%--------------------------------------------------------------------------
% remove cells outside volume
irem0 = any(finalpts(:,1:3)<1, 2);
iremx = finalpts(:,1) > size(annvol, 2);
iremy = finalpts(:,2) > size(annvol, 1);
iremz = finalpts(:,3) > size(annvol, 3);
irem  = irem0 | iremx | iremy | iremz;
finalpts(irem, :) = [];
igood(irem)       = [];
%--------------------------------------------------------------------------
% remove cells outside brain regions (root)
indcells = sub2ind(size(annvol), finalpts(:,2), finalpts(:,1), finalpts(:,3));
irem     = annvol(indcells) == 0;
finalpts(irem, :) = [];
igood(irem)       = [];
%--------------------------------------------------------------------------
finalpts = finalpts(igood, :);
fprintf('Retained %2.2f%% of detected cells in the brain\n', mean(igood)*100)
%--------------------------------------------------------------------------

% Nclust   = 20;
% Xuse     = finalpts(:,4:6);
% Xuse     = (Xuse - median(Xuse))./robustStd(Xuse);
% [aa, bb] = pca(Xuse);
% warning off;
% [idx, C] = kmeans(bb(:,1:3), Nclust, 'Replicates',10);
% warning on;
% %%
% sizethres = median(finalpts(:, 5)) + robustStd(finalpts(:, 5));
% elipsinds = accumarray(idx, finalpts(:, 6), [Nclust 1], @median);
% intinds   = accumarray(idx, finalpts(:, 4), [Nclust 1], @median);
% sizeinds  = accumarray(idx, finalpts(:, 5), [Nclust 1], @median);
% numinds   = accumarray(idx, 1, [Nclust 1], @sum);
% candinds  = intinds < median(finalpts(:, 4)) & sizeinds > sizethres;
% iremfin = ismember(idx, find(candinds));
% finalpts2  = finalpts;
% finalpts2(intremove, :) = [];
% % clf; scatter3(finalpts(:,1),finalpts(:,2),finalpts(:,3),1);hold on;
% clf;
% scatter3(finalpts(iremfin,1),finalpts(iremfin,2),finalpts(iremfin,3),1); axis equal; 
 %%
%--------------------------------------------------------------------------
end
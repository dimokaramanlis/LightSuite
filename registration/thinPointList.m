function keptIndices = thinPointList(cpts, thindist)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


% Pre-allocate memory for speed (worst case: keep all points)
Npoints        = size(cpts, 1);
keptIndices    = false(Npoints, 1); 
Dmat           = squareform(pdist(single(cpts)));
% We always keep the first point to start
keptIndices(1) = true;
%--------------------------------------------------------------------------
% We keep a separate list of kept points to vectorise the distance check
% (Initialize with the first point)

for ii = 2:Npoints
    
    distcurr = Dmat(ii, :);
    dists    = distcurr(keptIndices);
    
    % If the minimum distance to any kept point is greater than threshold
    if all(dists >= thindist)
        keptIndices(ii) = true;
    end
end
%--------------------------------------------------------------------------

end
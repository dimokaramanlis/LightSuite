function ptswarped = warpCordAtlas(pts, warppath)
%UNTITLED13 Summary of this function goes here
%   Detailed explanation goes here

ptsall = cell(numel(warppath), 1);
for inew = 1:numel(warppath)
    currids = pts(:, 2) == warppath(inew);
    ptsall{inew} = [pts(currids, 1) inew*ones(nnz(currids), 1) pts(currids, 3)];
end
ptswarped = cat(1, ptsall{:});
% pcshow(ptswarped)

end
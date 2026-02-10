function tvwarped = warpCordAtlasImage(tv, warppath)
%UNTITLED13 Summary of this function goes here
%   Detailed explanation goes here

idxwarp  = round(warppath);
idxwarp(idxwarp < 1) = 1;
idxwarp(idxwarp>size(tv, 1)) = size(tv, 1);
tvwarped = tv(idxwarp, :, :);

end
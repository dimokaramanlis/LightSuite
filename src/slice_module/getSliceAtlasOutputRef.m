function [Rout, yworld, ypix, nlow] = getSliceAtlasOutputRef(Rvolume, Ratlas, tform, nfac, dy)
% getSliceAtlasOutputRef returns the output view used to resample the
% rigidly transformed atlas into sample slice space.
%
% The legacy AP (Y) range covers the sample +- nfac slices. Here it is
% extended (keeping the legacy pixel grid) so that it also contains the
% whole transformed atlas. This lets the user pick any AP position in the
% control point GUI, even if the automated alignment put slices in the
% wrong AP location.
%
% nlow is the number of pixels prepended to the legacy range. Atlas slice
% indices stored in control points use the legacy convention, i.e.
% index_saved = index_in_Rout - nlow.

ylegacy   = [Rvolume.YWorldLimits(1) - dy*nfac, Rvolume.YWorldLimits(2) + dy*nfac];
ypixlegacy = ceil(range(ylegacy));
pxy       = range(ylegacy)/ypixlegacy;

[~, yatlas] = outputLimits(tform, Ratlas.XWorldLimits, Ratlas.YWorldLimits, Ratlas.ZWorldLimits);
nlow      = max(0, ceil((ylegacy(1) - yatlas(1))/pxy));
nhigh     = max(0, ceil((yatlas(2) - ylegacy(2))/pxy));

yworld    = [ylegacy(1) - nlow*pxy, ylegacy(2) + nhigh*pxy];
ypix      = ypixlegacy + nlow + nhigh;
Rout      = imref3d([ypix, Ratlas.ImageSize(2:3)], Rvolume.XWorldLimits, yworld, Rvolume.ZWorldLimits);
end

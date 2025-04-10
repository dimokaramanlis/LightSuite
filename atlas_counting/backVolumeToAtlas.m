function medianoverareas = backVolumeToAtlas(inputvol, trstruct)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%--------------------------------------------------------------------------
% at first, points should be resampled to the registration resolution

% we permute the volume to match atlas
volume  = permute(inputvol, trstruct.how_to_perm);
%--------------------------------------------------------------------------
volumereg = transformix(volume,trstruct.tform_bspline_samp20um_to_atlas_20um_px,...
    'movingscale', 0.02*[1 1 1]);
%--------------------------------------------------------------------------
volumereg          = uint16(abs(volumereg));
Rmoving            = imref3d(size(volumereg));
Rfixed             = imref3d(trstruct.atlassize);
registeredvolume   = imwarp(volumereg, Rmoving, trstruct.tform_affine_samp20um_to_atlas_10um_px,...
    'OutputView',Rfixed);
%--------------------------------------------------------------------------
% we reduce the signal to atlas areas
fprintf('Calculating background fluoresence in atlas coords... '); tic;

allen_atlas_path = fileparts(which('template_volume_10um.npy'));
av               = readNPY(fullfile(allen_atlas_path,'annotation_volume_10um_by_index.npy'));
% tv     = readNPY(fullfile(allen_atlas_path,'template_volume_10um.npy'));
Ngroups          = max(av, [], 'all');

Npxlr            = size(av,3)/2;
medianoverareas  = nan(Ngroups, 2, 'single');

for iside = 1:2
    istart = (iside - 1) * Npxlr + 1;
    iend   = istart + Npxlr - 1;

    sideav    = reshape(av(:, :, istart:iend), [], 1);
    sidevals  = reshape(registeredvolume(:, :, istart:iend), [], 1);
    ikeep     = sideav>1;
    medareas  = single(accumarray(sideav(ikeep), sidevals(ikeep), [Ngroups 1], @median));

    backlevel                 = single(median(sidevals(~ikeep)));
    medianoverareas(:, iside) = (medareas - backlevel)./backlevel;
end
fprintf('Done! Took %2.2f s\n', toc)

%--------------------------------------------------------------------------
end
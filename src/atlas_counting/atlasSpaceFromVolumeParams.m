function registeredvolume = atlasSpaceFromVolumeParams(inputvol, trstruct)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
%--------------------------------------------------------------------------
volume  = permuteBrainVolume(inputvol, trstruct.how_to_perm);
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
end
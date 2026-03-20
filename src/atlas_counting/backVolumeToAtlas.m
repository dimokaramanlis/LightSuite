function [medianoverareas, areaidx] = backVolumeToAtlas(inputvol, trstruct)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%--------------------------------------------------------------------------
registeredvolume = atlasSpaceFromVolumeParams(inputvol, trstruct);
%--------------------------------------------------------------------------
% we reduce the signal to atlas areas
fprintf('Calculating background fluoresence in atlas coords... '); tic;

allen_atlas_path = fileparts(which('annotation_10.nii.gz'));
av               = niftiread(fullfile(allen_atlas_path, 'annotation_10.nii.gz'));
parcelinfo       = readtable(fullfile(allen_atlas_path, 'parcellation_to_parcellation_term_membership.csv'));
areaidx          = unique(parcelinfo.parcellation_index);
Ngroups          = numel(areaidx);
Nforaccum        = max(av, [], 'all') + 1;

Npxlr            = size(av,3)/2;
medianoverareas  = nan(Ngroups, 2, 'single');

for iside = 1:2
    istart = (iside - 1) * Npxlr + 1;
    iend   = istart + Npxlr - 1;

    sideav    = reshape(av(:, :, istart:iend), [], 1);
    sidevals  = reshape(registeredvolume(:, :, istart:iend), [], 1);
    ikeep     = sideav>0;
    medareas  = single(accumarray(sideav(ikeep)+1, sidevals(ikeep), [Nforaccum 1], @median));
    medareas  = medareas(areaidx+1);

    backlevel                 = single(median(sidevals(~ikeep)));
    medareas(1)               = backlevel;
    medianoverareas(:, iside) = (medareas - backlevel)./backlevel;
end
fprintf('Done! Took %2.2f s\n', toc)

%--------------------------------------------------------------------------
end
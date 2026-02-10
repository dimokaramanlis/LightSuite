function [medianoverareas, areaidx] = backSliceVolumeToAtlas(inputvolpath, trstruct)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%--------------------------------------------------------------------------
fprintf('Loading data in memory... '); tic;
slicevol             = loadLargeSliceVolume(inputvolpath, 1);
fprintf('Done! Took %2.2f s\n', toc); 
%--------------------------------------------------------------------------
% we reduce the signal to atlas areas
fprintf('Calculating background fluoresence in atlas coords... '); tic;

allen_atlas_path = fileparts(which('annotation_10.nii.gz'));
av               = niftiread(fullfile(allen_atlas_path, 'annotation_10.nii.gz'));
av               = av(trstruct.atlasaplims(1):trstruct.atlasaplims(2), :, :);
parcelinfo       = readtable(fullfile(allen_atlas_path, 'parcellation_to_parcellation_term_membership.csv'));
areaidx          = unique(parcelinfo.parcellation_index);
Ngroups          = numel(areaidx);

% tv = niftiread(fullfile(allen_atlas_path,'average_template_10.nii.gz'));

Nforaccum        = max(av, [], 'all') + 1;

Npxlr            = size(av,3)/2;
medianoverareas  = nan(Ngroups, 2, 'single');

for iside = 1:2
    istart = (iside - 1) * Npxlr + 1;
    iend   = istart + Npxlr - 1;

    sideav    = reshape(av(:, :, istart:iend), [], 1);
    sidevals  = reshape(slicevol(:, :, istart:iend), [], 1);
    ikeep     = sideav>0;
    medareas  = single(accumarray(sideav(ikeep)+1, sidevals(ikeep), [Nforaccum 1], @median));
    medareas  = medareas(areaidx+1);

    backlevel                 = single(median(sidevals(~ikeep)));
    medareas(1)               = backlevel;
    if backlevel > 0
        medianoverareas(:, iside) = (medareas - backlevel)./backlevel;
    else
        medianoverareas(:, iside) = medareas;
    end
end
fprintf('Done! Took %2.2f s\n', toc)
%--------------------------------------------------------------------------
end
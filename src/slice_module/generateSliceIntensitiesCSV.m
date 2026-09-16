function allmedians = generateSliceIntensitiesCSV(sliceinfo, varargin)
%GENERATESLICEINTENSITIESCSV Per-region median intensity from a registered slice volume.
%
%   ALLMEDIANS = GENERATESLICEINTENSITIESCSV(SLICEINFO) loads the atlas-space
%   volume produced by generateRegisteredSliceVolume, overlays the Allen Brain
%   Atlas annotation, and computes per-region median intensities for every
%   channel.  Results are saved as chanNN_intensities.mat (always) and
%   chanNN_intensities.csv (optional) in <procpath>/volume_registered/.
%
%   Optional name-value pair:
%       'writetocsv'  - logical, default false
%
%   Output:
%       allmedians    - Ngroups x 2 x Nchans single array
%                       (dim 2: 1 = right hemisphere, 2 = left hemisphere)
%--------------------------------------------------------------------------
p = inputParser;
addRequired(p, 'sliceinfo');
addParameter(p, 'writetocsv', false, @(x) islogical(x) || isscalar(x));
parse(p, sliceinfo, varargin{:});
writetocsv = logical(p.Results.writetocsv);

atlasres_um = 10; % µm per voxel in the 10-µm Allen Atlas

registerpath = fullfile(sliceinfo.procpath, 'volume_registered');
if ~exist(registerpath, 'dir')
    error('Registered volume not found at %s.\nRun generateRegisteredSliceVolume first.', registerpath);
end
%--------------------------------------------------------------------------
% Load registered atlas-space volume: [NAP, NDV, Nchans, NLR]
fprintf('Loading registered volume... '); tic;
regvol = loadLargeSliceVolume(registerpath);
if ndims(regvol) == 3
    regvol = reshape(regvol, [size(regvol, 1), size(regvol, 2), 1, size(regvol, 3)]);
end
[~, ~, Nchans, ~] = size(regvol);
fprintf('Done! Took %2.2f s\n', toc);
%--------------------------------------------------------------------------
% Load regopts to get atlasaplims (AP crop used during registration)
ro = load(fullfile(sliceinfo.procpath, 'regopts.mat'));
%--------------------------------------------------------------------------
% Load atlas annotation at 10 µm and crop to AP range
fprintf('Loading atlas annotation... '); tic;
allen_atlas_path = fileparts(which('annotation_10.nii.gz'));
av = niftiread(fullfile(allen_atlas_path, 'annotation_10.nii.gz'));
av = av(ro.atlasaplims(1):ro.atlasaplims(2), :, :);

% Resize annotation to match registered volume's spatial dimensions [AP, DV, LR]
regvol_spatsize = size(regvol, [1 2 4]);
if ~isequal(size(av), regvol_spatsize)
    av = uint32(imresize3(single(av), regvol_spatsize, 'nearest'));
end
fprintf('Done! Took %2.2f s\n', toc);
%--------------------------------------------------------------------------
% Load parcellation info (same as generateRegisteredBrainVolumes)
allen_atlas_parcel_path = fileparts(which('parcellation_to_parcellation_term_membership.csv'));
parcelinfo  = readtable(fullfile(allen_atlas_parcel_path, ...
    'parcellation_to_parcellation_term_membership.csv'));

substridx        = strcmp(parcelinfo.parcellation_term_set_name, 'substructure');
[areaidx, ib]    = unique(parcelinfo.parcellation_index(substridx));
namessub         = parcelinfo.parcellation_term_name(substridx);
stridx           = strcmp(parcelinfo.parcellation_term_set_name, 'structure');
[~, ibstr]       = unique(parcelinfo.parcellation_index(stridx));
namesstruct      = parcelinfo.parcellation_term_name(stridx);
dividx           = strcmp(parcelinfo.parcellation_term_set_name, 'division');
[~, ibdiv]       = unique(parcelinfo.parcellation_index(dividx));
namesdiv         = parcelinfo.parcellation_term_name(dividx);

Ngroups   = numel(areaidx);
Nforaccum = double(max(av, [], 'all')) + 1;
Npxlr     = size(av, 3) / 2; % LR midline split (dim 3 = left-right axis)
%--------------------------------------------------------------------------
fprintf('Computing per-region intensities...\n'); proctic = tic;

allmedians = nan(Ngroups, 2, Nchans, 'single');

for ichan = 1:Nchans
    medianoverareas = nan(Ngroups, 2, 'single');
    volumeoverareas = nan(Ngroups, 2, 'single');

    vol_ichan = squeeze(regvol(:, :, ichan, :)); % [NAP, NDV, NLR]

    for iside = 1:2
        istart = (iside - 1) * Npxlr + 1;
        iend   = istart + Npxlr - 1;

        av_side  = reshape(double(av(:, :, istart:iend)), [], 1);
        vol_side = reshape(single(vol_ichan(:, :, istart:iend)), [], 1);
        ikeep    = av_side > 0;

        medareas = accumarray(av_side(ikeep) + 1, vol_side(ikeep), [Nforaccum 1], @median);
        medareas = medareas(areaidx + 1);

        backlevel              = single(median(vol_side(~ikeep)));
        medareas(areaidx == 0) = backlevel;
        medianoverareas(:, iside) = medareas;

        volareas = accumarray(av_side + 1, 1, [Nforaccum 1], @sum);
        volareas = volareas(areaidx + 1);
        volumeoverareas(:, iside) = volareas * (atlasres_um * 1e-3)^3;
    end

    allmedians(:, :, ichan) = medianoverareas;

    fmatname = fullfile(registerpath, sprintf('chan%02d_intensities.mat', ichan));
    save(fmatname, 'medianoverareas', 'areaidx', 'volumeoverareas');

    if writetocsv
        currtable = array2table([areaidx, medianoverareas, volumeoverareas], ...
            'VariableNames', {'parcellation_index', ...
            'RightSideIntensity', 'LeftSideIntensity', ...
            'RightSideVolume[mm3]', 'LeftSideVolume[mm3]'});
        currtable = addvars(currtable, namessub(ib), namesstruct(ibstr), namesdiv(ibdiv), ...
            'NewVariableNames', {'name','structure','division'}, 'Before', 'parcellation_index');
        fsavename = fullfile(registerpath, sprintf('chan%02d_intensities.csv', ichan));
        writetable(currtable, fsavename, 'Delimiter', ';');
    end

    fprintf('Channel %d/%d done. Time %2.2f s.\n', ichan, Nchans, toc(proctic));
end
%--------------------------------------------------------------------------
end

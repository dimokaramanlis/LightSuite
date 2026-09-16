function sliceCellCountsToCSV(procpath, varargin)
%SLICECELLCOUNTSTOCSV Count atlas-space slice cells per brain region and save as mat/CSV.
%
%   SLICECELLCOUNTSTOCSV(PROCPATH) loads per-channel atlas cell coordinate
%   files produced by slicePointsToAtlas, counts cells in each Allen Brain
%   Atlas region (both hemispheres), and saves chanNN_cellcounts.mat and
%   (optionally) chanNN_cellcounts.csv to <procpath>/volume_registered/.
%
%   File search order in PROCPATH:
%       1. chan*_cell_locations_atlas.mat  (per-channel, from multi-channel workflow)
%       2. cell_locations_atlas.mat        (legacy single-channel file)
%
%   Optional name-value pair:
%       'writetocsv'  - logical, default true
%--------------------------------------------------------------------------
p = inputParser;
addRequired(p, 'procpath', @(x) ischar(x) || isstring(x));
addParameter(p, 'writetocsv', true, @(x) islogical(x) || isscalar(x));
parse(p, procpath, varargin{:});
writetocsv = logical(p.Results.writetocsv);

registerpath = fullfile(procpath, 'volume_registered');
makeNewDir(registerpath);
%--------------------------------------------------------------------------
% Load Allen Atlas annotation (full, uncropped — atlas coords from
% slicePointsToAtlas already account for atlasaplims offset)
fprintf('Loading Allen Atlas annotation...\n');
allen_atlas_path = fileparts(which('annotation_10.nii.gz'));
av = niftiread(fullfile(allen_atlas_path, 'annotation_10.nii.gz'));
%--------------------------------------------------------------------------
% Load parcellation info (same pattern as transformPointsToAtlas)
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
%--------------------------------------------------------------------------
% Find atlas cell location files
files = dir(fullfile(procpath, 'chan*_cell_locations_atlas.mat'));
if isempty(files)
    files = dir(fullfile(procpath, 'cell_locations_atlas.mat')); % legacy fallback
end
if isempty(files)
    fprintf('No cell_locations_atlas.mat files found in %s.\n', procpath);
    return;
end
%--------------------------------------------------------------------------
for i = 1:numel(files)
    currfile = fullfile(files(i).folder, files(i).name);
    fprintf('Processing %s... ', files(i).name); savetic = tic;

    dat    = load(currfile);
    fnames = fieldnames(dat);
    atlasptcoords = dat.(fnames{1}); % first field, regardless of variable name

    % Sanitize and remove out-of-brain points
    [cleanpts, ~] = sanitizeCellCoords(atlasptcoords, av);

    % Count cells per leaf region, both hemispheres
    [areacounts, areavols, ~] = groupCellsIntoLeafRegions(cleanpts, av, areaidx);

    % Extract channel number from filename (e.g., "chan02_cell_locations_atlas.mat")
    tokens = regexp(files(i).name, 'chan(\d+)', 'tokens');
    if isempty(tokens)
        ichan = i;
    else
        ichan = str2double(tokens{1}{1});
    end

    % Save .mat
    fmatname = fullfile(registerpath, sprintf('chan%02d_cellcounts.mat', ichan));
    save(fmatname, 'areacounts', 'areaidx', 'areavols');

    % Save .csv
    if writetocsv
        currtable = array2table([areaidx, areacounts, areavols], ...
            'VariableNames', ...
            {'parcellation_index','LeftSideCount','RightSideCount','TotalVolume[mm3]'});
        currtable = addvars(currtable, namessub(ib), namesstruct(ibstr), namesdiv(ibdiv), ...
            'NewVariableNames', {'name','structure','division'}, 'Before', 'parcellation_index');
        fsavename = fullfile(registerpath, sprintf('chan%02d_cellcounts.csv', ichan));
        writetable(currtable, fsavename, 'Delimiter', ';');
    end

    fprintf('Done! Took %2.2f s.\n', toc(savetic));
end
%--------------------------------------------------------------------------
end

function varargout = transformPointsToAtlas(input_data, varargin)
%TRANSFORMPOINTSTOATLAS Transforms cell coordinates to atlas space and counts them.
%
%   FINALPTS = TRANSFORMPOINTSTOATLAS(INPUT_PTS, 'transform_params', TRSTRUCT)
%   transforms the N x M array of INPUT_PTS directly. N is the number of
%   detected points, and M must be >= 3. The first 3 columns represent the
%   [x, y, z] spatial coordinates in the original sample space. Any additional
%   columns (e.g., intensity, equivalent diameter, elongation) are purely
%   descriptive and will be carried over to the output unmodified. TRSTRUCT
%   must contain the transformation parameters.
%
%   TRANSFORMPOINTSTOATLAS(SAVEPATH) searches the directory SAVEPATH for point
%   files, transforms them, saves the result as '*_atlas.mat' inside the
%   'volume_registered' directory, and calculates cell counts per brain region
%   using the Allen Brain Atlas.
%
%   TRANSFORMPOINTSTOATLAS(FILEPATH, 'savepath', SAVEPATH) does the same for one
%   named point file that lives outside the registration folder.
%
%   Point file formats
%   ------------------
%     .mat  LightSuite detections: a 'cell_locations' array [N x M], M >= 3,
%           columns [x y z, descriptors...]. If the file also holds
%           'cell_images' (written when opts.savecellimages is on) the CNN
%           artifact classifier can be run on it, see 'network' below.
%     .csv  the same array as text, i.e. what writematrix(cell_locations, ...)
%           writes. A header line is allowed and skipped.
%     .xml  an ImageJ / Fiji "Cell Counter" marker file. Every <Marker_Type>
%           block becomes one point set of [x y z] columns.
%
%   In folder mode the search patterns are '*cell_locations_sample.mat',
%   '*cell_locations_sample.csv' and '*.xml'.
%
%   Channels
%   --------
%   Outputs are written per channel, so every point set has to name the channel
%   it came from. It is taken, in this order, from the 'channel' option, from a
%   'chan_3_' / 'chan03_' / 'channel3_' tag in the file name, or - for XML only -
%   from the <Type> of the marker block. A file that names no channel at all
%   falls back to its position in the list, as it always has, but warns that it
%   guessed. The resolved mapping is printed for every point set either way, and
%   two point sets claiming the same channel warn too, since their per-channel
%   statistics would overwrite each other.
%
%   Marker types are the Cell Counter plugin's categories, NOT image channels;
%   nothing in the file records the channel. Map them yourself when they differ:
%
%       transformPointsToAtlas('D:\data\cells.xml', 'savepath', savepath, ...
%           'channel', [2 3]);   % marker Type 1 -> chan 2, Type 2 -> chan 3
%
%   Cell classification
%   -------------------
%   TRANSFORMPOINTSTOATLAS(..., 'network', NET) runs the CNN artifact classifier
%   before the transform, so only candidates it accepts as real cells reach atlas
%   space and the regional counts. NET is a path to a .mat saved by
%   demos/trainClassificationNetwork.m, or the network object itself.
%
%   The classification of a file is computed once and cached next to it as
%   '<name>_classification.mat'. A later run reuses that file instead of
%   re-running the network, unless 'reclassify' is true or the cached result no
%   longer matches the number of detections.
%
%   Classification needs the detector's cell images, so it only applies to .mat
%   point files saved with opts.savecellimages = true; CSV and XML point sets are
%   transformed unfiltered, with a warning.
%
%   Inputs:
%       input_data           - (char/string) Directory path, or path to a single
%                              point file, OR
%                              (numeric) N x M array of points (N points,
%                              M dimensions where cols 1:3 are x, y, z
%                              in sample space, and cols 4:end are single
%                              number descriptors like intensity/diameter).
%
%   Optional Name-Value Parameters:
%       transform_params     - (struct) Transformation structure. Required
%                              if input_data is a numeric array.
%       savepath             - (char/string) Folder holding transform_params.mat
%                              and regopts.mat. Required when input_data is a
%                              single file outside that folder; defaults to the
%                              file's own folder.
%       writetocsv           - (logical) If true, writes data to CSV.
%       channel              - (numeric) Channel of each point set. A scalar
%                              applies to all of them; a vector must have one
%                              entry per point set of the file being read.
%       network              - Trained cell classifier (path or object). Empty
%                              (default) skips classification.
%       reclassify           - (logical) Re-run the network even when a cached
%                              classification exists. Default false.
%       coordoffset          - (1x3) Offset added to XML markers, see
%                              readCellCounterXML. Default [1 1 0].
%
%   Outputs:
%       finalpts             - (single/double) Array of transformed points, or a
%                              cell array of them when several point sets were
%                              processed. Returned if an output is requested.
%
%   See also READCELLLOCATIONSFILE, READCELLCOUNTERXML, CLASSIFYDETECTEDCELLS,
%   GENERATEREGISTEREDBRAINVOLUMES.
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% 1. Parse Inputs and establish defaults
%--------------------------------------------------------------------------
p = inputParser;
addRequired(p, 'input_data', @(x) ischar(x) || isstring(x) || isnumeric(x));
addParameter(p, 'transform_params', [], @isstruct);
addParameter(p, 'writetocsv', [], @(x) islogical(x) || isscalar(x));
addParameter(p, 'savepath', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'channel', [], @(x) isnumeric(x) && all(x == round(x)) && all(x >= 0));
addParameter(p, 'network', []);
addParameter(p, 'reclassify', false, @(x) islogical(x) || isscalar(x));
addParameter(p, 'coordoffset', [1 1 0], @(x) isnumeric(x) && numel(x) == 3);

parse(p, input_data, varargin{:});

if isnumeric(input_data)
    % Points mode (direct transformation without statistics)
    trstruct = p.Results.transform_params;
    if isempty(trstruct)
        error('''transform_params'' structure must be provided when passing points directly.');
    end
    finalpts = coreTransform(input_data, trstruct);
    if nargout > 0
        varargout{1} = finalpts;
    end
    return;
end

%--------------------------------------------------------------------------
% 2. Resolve the registration folder and the list of point files
%--------------------------------------------------------------------------
inpath = char(input_data);

if isfolder(inpath)
    savepath  = inpath;
    pointpaths = collectPointFiles(savepath);
    if isempty(pointpaths)
        fprintf('No valid cell location files found to process.\n');
        if nargout > 0; varargout{1} = []; end
        return;
    end
elseif isfile(inpath)
    savepath = char(p.Results.savepath);
    if isempty(savepath)
        savepath = fileparts(inpath);
    end
    pointpaths = {inpath};
else
    error('transformPointsToAtlas:badInput', ...
        '%s is neither a folder nor a file.', inpath);
end

fprintf('Running in folder mode. Target: %s\n', savepath);

trfile  = fullfile(savepath, 'transform_params.mat');
optfile = fullfile(savepath, 'regopts.mat');
if ~exist(trfile, 'file') || ~exist(optfile, 'file')
    error(['transform_params.mat or regopts.mat not found in %s. Pass the ' ...
        'registration folder with ''savepath'' if the point file lives ' ...
        'elsewhere.'], savepath);
end

trstruct   = load(trfile);
loadedOpts = load(optfile);
opts       = loadedOpts.opts;

% Determine writetocsv fallback
writetocsv = false;
if ~isempty(p.Results.writetocsv)
    writetocsv = p.Results.writetocsv;
elseif isfield(opts, 'writetocsv')
    writetocsv = opts.writetocsv;
end

% Load the classifier once, if one was asked for
netin      = p.Results.network;
useclassif = ~isempty(netin);
netinfo    = struct('name', '', 'path', '', 'accuracy', NaN);
if useclassif
    [net, netinfo] = loadCellClassifierNet(netin);
    fprintf('Cell classifier: %s\n', netinfo.name);
end

% Setup registration output directory
registerpath = fullfile(savepath, 'volume_registered');
if ~exist(registerpath, 'dir')
    mkdir(registerpath);
end

% Load Allen Atlas and Parcellation info once
fprintf('Loading Allen Atlas and parcellation data...\n');
atlasData = loadAtlasData();
groupinds = atlasData.areaidx;

%--------------------------------------------------------------------------
% 3. Read every point file into channel-tagged datasets
%--------------------------------------------------------------------------
readopts = struct('channel', p.Results.channel, ...
                  'coordoffset', p.Results.coordoffset, 'verbose', true);

datasets = [];
for i = 1:numel(pointpaths)
    fprintf('Reading %s...\n', pointpaths{i});
    % keep the historical behaviour for files that name no channel: fall back
    % to their position in the list, but say so loudly
    readopts.defaultchannel = i;
    try
        ds = readCellLocationsFile(pointpaths{i}, readopts);
    catch ME
        [~, ~, thisext] = fileparts(pointpaths{i});
        if strcmpi(thisext, '.xml')
            % folder mode globs every .xml; an unrelated one is not an error
            warning('transformPointsToAtlas:skipXml', ...
                'Skipping %s: %s', pointpaths{i}, ME.message);
            continue;
        end
        rethrow(ME);
    end
    if isempty(datasets)
        datasets = ds;
    else
        datasets = [datasets, ds]; %#ok<AGROW>
    end
end

if isempty(datasets)
    fprintf('No valid cell location files found to process.\n');
    if nargout > 0; varargout{1} = []; end
    return;
end

checkChannelClashes(datasets);

%--------------------------------------------------------------------------
% 4. Process each dataset
%--------------------------------------------------------------------------
all_final_pts = cell(numel(datasets), 1);
for i = 1:numel(datasets)
    ds    = datasets(i);
    ichan = ds.channel;
    fprintf('Processing %s (channel %d)... \n', ds.label, ichan); savetic = tic;

    inputpts = ds.points;
    nraw     = size(inputpts, 1);
    classres = [];

    % 1. Optional CNN classification: keep only the candidates it calls cells
    if useclassif
        [inputpts, classres] = applyCellClassifierCached(inputpts, ds, net, netinfo, ...
            savepath, p.Results.reclassify);
    end

    if isempty(inputpts)
        warning('transformPointsToAtlas:noPointsLeft', ...
            'No points left for %s after classification; skipping.', ds.label);
        continue;
    end

    % 2. Transform Points
    atlasptcoords = coreTransform(inputpts, trstruct);

    % 3. Sanitize and remove invalid points
    [cleanatlaspts, badpts] = sanitizeCellCoords(atlasptcoords, atlasData.av);

    % 4. Group into brain regions and calculate stats
    [areacounts, areavols, catids] = groupCellsIntoLeafRegions(...
        cleanatlaspts, atlasData.av, groupinds);

    % 5. Save locations in the volume_registered folder
    outname   = sprintf('%s_atlas.mat', ds.outbase);
    fsavename = fullfile(registerpath, outname);
    sourcefile = ds.sourcefile;
    if isempty(classres)
        save(fsavename, 'atlasptcoords', 'cleanatlaspts', 'badpts', 'catids', ...
            'ichan', 'sourcefile');
    else
        classification = rmfield(classres, {'labels', 'scores'});
        save(fsavename, 'atlasptcoords', 'cleanatlaspts', 'badpts', 'catids', ...
            'ichan', 'sourcefile', 'classification');
    end

    if writetocsv
        loc_csvname = strrep(outname, '.mat', '.csv');
        writematrix(atlasptcoords, fullfile(registerpath, loc_csvname));
    end

    % 6. Save statistics
    saveCellStats(registerpath, ichan, areacounts, areavols, atlasData, writetocsv);

    all_final_pts{i} = atlasptcoords;
    if isempty(classres)
        fprintf('  -> Done! Channel %d, %d points, completed in %2.2f s.\n', ...
            ichan, size(atlasptcoords, 1), toc(savetic));
    else
        fprintf('  -> Done! Channel %d, %d/%d points kept by the classifier, completed in %2.2f s.\n', ...
            ichan, size(inputpts, 1), nraw, toc(savetic));
    end
end

if nargout > 0
    if isscalar(all_final_pts)
        varargout{1} = all_final_pts{1};
    else
        varargout{1} = all_final_pts;
    end
end
fprintf('All points successfully registered and tabulated!\n');

end

%==========================================================================
% Local Helper Functions
%==========================================================================

function pointpaths = collectPointFiles(savepath)
%COLLECTPOINTFILES Point files of every supported format inside savepath.
patterns   = {'*cell_locations_sample.mat', '*cell_locations_sample.csv', '*.xml'};
pointpaths = {};
for k = 1:numel(patterns)
    f = dir(fullfile(savepath, patterns{k}));
    for j = 1:numel(f)
        pointpaths{end+1} = fullfile(f(j).folder, f(j).name); %#ok<AGROW>
    end
end
pointpaths = unique(pointpaths, 'stable');
end

%--------------------------------------------------------------------------
function checkChannelClashes(datasets)
%CHECKCHANNELCLASHES Warn when two point sets would overwrite each other.
%   Per-channel statistics are written as chanNN_cellcounts.mat, so two point
%   sets claiming the same channel silently overwrite one another.
chans = [datasets.channel];
[uc, ~, ic] = unique(chans);
for k = 1:numel(uc)
    hits = find(ic == k);
    if numel(hits) > 1
        warning('transformPointsToAtlas:duplicateChannel', ...
            ['Channel %d is claimed by %d point sets (%s). Their per-channel ' ...
             'statistics will overwrite each other - pass ''channel'' to map ' ...
             'them onto distinct channels.'], ...
            uc(k), numel(hits), strjoin({datasets(hits).label}, ', '));
    end
end
end


%--------------------------------------------------------------------------
function finalpts = coreTransform(inputpts, trstruct)
    regsize_mm = trstruct.atlasres * 2 * 1e-3;

    pts = (inputpts(:, 1:3) - 1) .* trstruct.ori_pxsize * 1e-3;
    pts = pts(:, [2 1 3]);

    if isfield(trstruct, 'ori_size')
        phys_size_yxz = (trstruct.ori_size - 1) .* trstruct.ori_pxsize([2 1 3]) * 1e-3;
    else
        error('trstruct must contain ''ori_size'' for axis flips.');
    end

    perm_order   = abs(trstruct.how_to_perm);
    pts_permuted = zeros(size(pts), 'like', pts);

    for dim = 1:3
        orig_dim = perm_order(dim);
        pts_permuted(:, dim) = pts(:, orig_dim);
        if trstruct.how_to_perm(dim) < 0
            dim_max_size = phys_size_yxz(orig_dim);
            pts_permuted(:, dim) = dim_max_size - pts_permuted(:, dim);
        end
    end

    pts = pts_permuted(:, [2 1 3]) / regsize_mm;

    Dfield = transformix([], trstruct.tform_bspline_samp20um_to_atlas_20um_px);
    Dfield = permute(Dfield, [2 3 4 1]) / regsize_mm;
    [Sx, Sy, Sz, ~] = size(Dfield);

    Xgv = 1:Sx; Ygv = 1:Sy; Zgv = 1:Sz;
    dx = interpn(Xgv, Ygv, Zgv, Dfield(:,:,:,1), pts(:,1), pts(:,2), pts(:,3), 'linear');
    dy = interpn(Xgv, Ygv, Zgv, Dfield(:,:,:,2), pts(:,1), pts(:,2), pts(:,3), 'linear');
    dz = interpn(Xgv, Ygv, Zgv, Dfield(:,:,:,3), pts(:,1), pts(:,2), pts(:,3), 'linear');

    interpolated_displacements = -[dx, dy, dz];
    interpolated_displacements(isnan(interpolated_displacements)) = 0;

    finalpts = trstruct.tform_affine_samp20um_to_atlas_10um_px.transformPointsForward(pts + interpolated_displacements);
    finalpts = [finalpts, inputpts(:, 4:end)];
end

function atlasData = loadAtlasData()
    % Loads the annotation volume and parcellation CSV info
    allen_atlas_path        = fileparts(which('annotation_10.nii.gz'));
    av                      = niftiread(fullfile(allen_atlas_path, 'annotation_10.nii.gz'));
    allen_atlas_parcel_path = fileparts(which('parcellation_to_parcellation_term_membership.csv'));

    parcelinfo       = readtable(fullfile(allen_atlas_parcel_path, 'parcellation_to_parcellation_term_membership.csv'));

    substridx        = strcmp(parcelinfo.parcellation_term_set_name, 'substructure');
    [areaidx, ib]    = unique(parcelinfo.parcellation_index(substridx));
    namessub         = parcelinfo.parcellation_term_name(substridx);

    stridx           = strcmp(parcelinfo.parcellation_term_set_name, 'structure');
    [~, ibstr]       = unique(parcelinfo.parcellation_index(stridx));
    namesstruct      = parcelinfo.parcellation_term_name(stridx);

    dividx           = strcmp(parcelinfo.parcellation_term_set_name, 'division');
    [~, ibdiv]       = unique(parcelinfo.parcellation_index(dividx));
    namesdiv         = parcelinfo.parcellation_term_name(dividx);

    % Pack into a clean struct
    atlasData.av          = av;
    atlasData.parcelinfo  = parcelinfo;
    atlasData.areaidx     = areaidx;
    atlasData.namessub    = namessub(ib);
    atlasData.namesstruct = namesstruct(ibstr);
    atlasData.namesdiv    = namesdiv(ibdiv);
end

function saveCellStats(registerpath, ichan, areacounts, areavols, atlasData, writetocsv)
    % Saves the cell counts mirroring generateRegisteredBrainVolumes
    fmatname = fullfile(registerpath, sprintf('chan%02d_cellcounts.mat', ichan));
    areaidx  = atlasData.areaidx;
    save(fmatname, 'areacounts', 'areaidx', 'areavols');

    if writetocsv
        currtable = array2table([areaidx, areacounts, areavols], ...
            'VariableNames', ...
            {'parcellation_index', 'RightSideCount', 'LeftSideCount', 'TotalVolume[mm3]'});

        currtable = addvars(currtable, atlasData.namessub, atlasData.namesstruct, atlasData.namesdiv, ...
            'NewVariableNames', {'name', 'structure', 'division'}, 'Before', 'parcellation_index');

        fsavename = fullfile(registerpath, sprintf('chan%02d_cellcounts.csv', ichan));
        writetable(currtable, fsavename);
    end
end

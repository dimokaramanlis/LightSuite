function varargout = transformCordPointsToAtlas(input_data, varargin)
%TRANSFORMCORDPOINTSTOATLAS Move detected points into spinal cord atlas space.
%
%   FINALPTS = TRANSFORMCORDPOINTSTOATLAS(INPUT_PTS, 'transform_params', TRSTRUCT)
%   transforms the N x M array INPUT_PTS directly. Columns 1:3 are [x y z] in the
%   original sample, the coordinates cell detection produces; columns 4 and
%   beyond are descriptors and are carried over unchanged.
%
%   TRANSFORMCORDPOINTSTOATLAS(SAVEPATH) does the same for every point file in
%   the registration folder SAVEPATH, saves each result as '*_atlas.mat' inside
%   'volume_registered', and counts the points per atlas region and per spinal
%   segment. It is the cord counterpart of TRANSFORMPOINTSTOATLAS and reads the
%   same file formats through READCELLLOCATIONSFILE: LightSuite .mat detections,
%   the same array as .csv, and ImageJ "Cell Counter" .xml marker files, with the
%   same rules for working out which channel a point set belongs to.
%
%   TRANSFORMCORDPOINTSTOATLAS(FILEPATH, 'savepath', SAVEPATH) handles one named
%   point file that lives outside the registration folder.
%
%   What makes the cord different is the transform itself, not the bookkeeping:
%   the sample was permuted, cropped at both ends, and straightened slice by
%   slice before any atlas was fitted to it, so all of that has to be undone in
%   order. CORDPOINTSTOATLAS does that and documents the chain; counts come out
%   per region *and* per segment, since a cord region runs the whole length of
%   the cord and a single number for it would hide the rostrocaudal axis.
%
%   Optional name-value arguments
%     transform_params - transform struct; required when points are passed
%                        directly, otherwise read from SAVEPATH.
%     savepath         - folder holding transform_params.mat and regopts.mat.
%     writetocsv       - also write the counts as CSV (default from
%                        opts.writetocsv, else false).
%     channel          - channel of each point set, see READCELLLOCATIONSFILE.
%     network          - trained cell classifier (path or object) applied before
%                        the transform. Empty (default) skips classification.
%     reclassify       - re-run the network even when a cached result exists.
%     coordoffset      - 1x3 offset added to XML markers, default [1 1 0].
%
%   See also CORDPOINTSTOATLAS, GENERATEREGISTEREDCORDVOLUME,
%   TRANSFORMPOINTSTOATLAS, READCELLLOCATIONSFILE, GROUPCORDPOINTSINTOAREAS.

%--------------------------------------------------------------------------
% 1. Parse inputs
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
    trstruct = p.Results.transform_params;
    if isempty(trstruct)
        error('transformCordPointsToAtlas:noTransform', ...
            '''transform_params'' must be provided when passing points directly.');
    end
    finalpts = cordPointsToAtlas(input_data, trstruct);
    if nargout > 0
        varargout{1} = finalpts;
    end
    return;
end

%--------------------------------------------------------------------------
% 2. Resolve the registration folder and the point files
%--------------------------------------------------------------------------
inpath = char(input_data);

if isfolder(inpath)
    savepath   = inpath;
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
    error('transformCordPointsToAtlas:badInput', ...
        '%s is neither a folder nor a file.', inpath);
end

fprintf('Running in folder mode. Target: %s\n', savepath);

trfile = fullfile(savepath, 'transform_params.mat');
if ~exist(trfile, 'file')
    error('transformCordPointsToAtlas:noTransformFile', ...
        ['transform_params.mat not found in %s. Pass the registration folder ' ...
         'with ''savepath'' if the point file lives elsewhere.'], savepath);
end
trstruct = load(trfile);
opts     = loadRegOpts(savepath);

writetocsv = getOr(opts, 'writetocsv', false);
if ~isempty(p.Results.writetocsv)
    writetocsv = p.Results.writetocsv;
end

% classifier, loaded once
netin      = p.Results.network;
useclassif = ~isempty(netin);
netinfo    = struct('name', '', 'path', '', 'accuracy', NaN);
net        = [];
if useclassif
    [net, netinfo] = loadCellClassifierNet(netin);
    fprintf('Cell classifier: %s\n', netinfo.name);
end

registerpath = fullfile(savepath, 'volume_registered');
makeNewDir(registerpath);

%--------------------------------------------------------------------------
% 3. Atlas and its region/segment lookups, loaded once
%--------------------------------------------------------------------------
fprintf('Loading the spinal cord atlas...\n');
[~, av, parcelinfo, segmentinfo, atlasres] = loadSpinalCordAtlas();
grouping = cordAtlasGrouping(av, segmentinfo);

% names of every region at all three levels, so the counts identify a region the
% way the brain ones do rather than by atlas id alone
areahierarchy = cordAreaHierarchy(parcelinfo, grouping.avinds);

% the displacement field is the expensive part of the transform, so read it
% once and hand it to every point set
Dfield = transformix([], trstruct.tform_bspline_samp20um_to_atlas_20um_px);
Dfield = permute(Dfield, [2 3 4 1]) / (trstruct.registrationres(1) * 1e-3);

%--------------------------------------------------------------------------
% 4. Read every point file into channel-tagged datasets
%--------------------------------------------------------------------------
readopts = struct('channel', p.Results.channel, ...
                  'coordoffset', p.Results.coordoffset, 'verbose', true);

datasets = [];
for ii = 1:numel(pointpaths)
    fprintf('Reading %s...\n', pointpaths{ii});
    readopts.defaultchannel = ii;
    try
        ds = readCellLocationsFile(pointpaths{ii}, readopts);
    catch ME
        [~, ~, thisext] = fileparts(pointpaths{ii});
        if strcmpi(thisext, '.xml')
            % folder mode globs every .xml; an unrelated one is not an error
            warning('transformCordPointsToAtlas:skipXml', ...
                'Skipping %s: %s', pointpaths{ii}, ME.message);
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
% 5. Transform and count each dataset
%--------------------------------------------------------------------------
all_final_pts = cell(numel(datasets), 1);
for ii = 1:numel(datasets)
    ds    = datasets(ii);
    ichan = ds.channel;
    fprintf('Processing %s (channel %d)... \n', ds.label, ichan); savetic = tic;

    inputpts = ds.points;
    nraw     = size(inputpts, 1);
    classres = [];

    if useclassif
        [inputpts, classres] = applyCellClassifierCached(inputpts, ds, net, netinfo, ...
            savepath, p.Results.reclassify);
    end

    if isempty(inputpts)
        warning('transformCordPointsToAtlas:noPointsLeft', ...
            'No points left for %s after classification; skipping.', ds.label);
        continue;
    end

    atlasptcoords = cordPointsToAtlas(inputpts, trstruct, 'Dfield', Dfield);

    [cleanatlaspts, badpts] = sanitizeCellCoords(atlasptcoords, av);

    [areacounts, areavols, ptareas, ptsegments] = groupCordPointsIntoAreas(...
        cleanatlaspts, av, grouping, atlasres);

    outname    = sprintf('%s_atlas.mat', ds.outbase);
    fsavename  = fullfile(registerpath, outname);
    sourcefile = ds.sourcefile;
    if isempty(classres)
        save(fsavename, 'atlasptcoords', 'cleanatlaspts', 'badpts', ...
            'ptareas', 'ptsegments', 'ichan', 'sourcefile');
    else
        classification = rmfield(classres, {'labels', 'scores'});
        save(fsavename, 'atlasptcoords', 'cleanatlaspts', 'badpts', ...
            'ptareas', 'ptsegments', 'ichan', 'sourcefile', 'classification');
    end

    if writetocsv
        writematrix(atlasptcoords, fullfile(registerpath, strrep(outname, '.mat', '.csv')), 'Delimiter', ';');
    end

    saveCordCellStats(registerpath, ichan, areacounts, areavols, grouping, ...
        areahierarchy, writetocsv);

    all_final_pts{ii} = atlasptcoords;
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
% Local helpers
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
chans = [datasets.channel];
[uc, ~, ic] = unique(chans);
for k = 1:numel(uc)
    hits = find(ic == k);
    if numel(hits) > 1
        warning('transformCordPointsToAtlas:duplicateChannel', ...
            ['Channel %d is claimed by %d point sets (%s). Their per-channel ' ...
             'statistics will overwrite each other - pass ''channel'' to map ' ...
             'them onto distinct channels.'], ...
            uc(k), numel(hits), strjoin({datasets(hits).label}, ', '));
    end
end
end

%--------------------------------------------------------------------------
function saveCordCellStats(registerpath, ichan, areacounts, areavols, grouping, ...
    areahierarchy, writetocsv)
%SAVECORDCELLSTATS Per-region, per-segment counts, mirroring the brain outputs.

areaidx     = grouping.avinds;
segmentname = grouping.segnames;

fmatname = fullfile(registerpath, sprintf('chan%02d_cellcounts.mat', ichan));
save(fmatname, 'areacounts', 'areaidx', 'areavols', 'segmentname', 'areahierarchy');

if writetocsv
    currtable = cordAreaTable(areahierarchy, segmentname, ...
        struct('Count', areacounts, 'Volume_mm3', areavols));
    writetable(currtable, fullfile(registerpath, sprintf('chan%02d_cellcounts.csv', ichan)), 'Delimiter', ';');
end

end

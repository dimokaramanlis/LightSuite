function datasets = readCellLocationsFile(filepath, opts)
%READCELLLOCATIONSFILE Load detected cell positions from .mat, .csv or .xml.
%   DATASETS = READCELLLOCATIONSFILE(FILEPATH) reads one file of detected cell
%   positions and returns a struct array with one entry per point set in it,
%   each tagged with the image CHANNEL it belongs to. It is the input layer of
%   transformPointsToAtlas, so every format below reaches the transform in the
%   same shape.
%
%   Supported formats
%   -----------------
%   .mat  A LightSuite detection file holding a 'cell_locations' array
%         [N x M], M >= 3, columns [x y z, descriptors...] in sample space.
%         If the file also holds 'cell_images' (written when
%         opts.savecellimages is on) they are carried along, so the CNN
%         artifact classifier can be run on them later.
%
%   .csv  The same array written as text, i.e. what
%         writematrix(cell_locations, ...) produces. A header line is allowed
%         and skipped; columns must still be [x y z, descriptors...]. CSV
%         files carry no cell images, so they cannot be classified.
%
%   .xml  An ImageJ / Fiji "Cell Counter" marker file (readCellCounterXML).
%         Each <Marker_Type> block becomes one dataset with columns [x y z].
%         XML files carry no cell images either.
%
%   Channel assignment
%   ------------------
%   Every dataset must name the channel it came from, because the atlas-space
%   outputs are written per channel (chanNN_cellcounts.mat, ...). The channel
%   is resolved in this order:
%
%     1. opts.channel, if given (see below);
%     2. the file name: 'chan_3_...', 'chan03_...', 'channel3_...';
%     3. for XML only, the <Marker_Type> <Type> value;
%     4. opts.defaultchannel, if the caller supplied a fallback - used with a
%        warning, since it is a guess rather than something the data says.
%
%   With no fallback the function errors rather than guessing. Either way the
%   channel of every point set, and where it came from, is printed.
%
%   Input:
%     filepath   full path to a .mat, .csv or .xml file.
%     opts       (optional) struct:
%        .channel        channel number(s) to assign. A scalar applies to every
%                        dataset in the file; a vector must have one entry per
%                        dataset (for XML, one per <Marker_Type>, in file order).
%        .defaultchannel fallback channel when nothing else identifies one.
%                        Empty (default) makes an unidentifiable channel an
%                        error instead.
%        .coordoffset    1x3 offset for XML markers (see readCellCounterXML);
%                        default [1 1 0].
%        .verbose        default true.
%
%   Output:
%     datasets   struct array with fields
%        .points       [N x M] single, columns [x y z, descriptors...]
%        .channel      resolved channel number
%        .outbase      base name for the atlas-space output files
%        .sourcefile   filepath
%        .sourcetype   'mat' | 'csv' | 'xml'
%        .label        short human-readable description (for messages)
%        .cell_images  [N x nsigma] or [] if the format carries none
%        .imwindow     1x3 half-window of cell_images, or []
%        .markertype   XML <Type> value, or NaN
%
%   See also TRANSFORMPOINTSTOATLAS, READCELLCOUNTERXML.

if nargin < 2, opts = struct(); end
chanopt     = getOr(opts, 'channel',        []);
defaultchan = getOr(opts, 'defaultchannel', []);
coordoffset = getOr(opts, 'coordoffset',    [1 1 0]);
verbose     = getOr(opts, 'verbose',        true);

if ~isfile(filepath)
    error('readCellLocationsFile:fileNotFound', 'File not found: %s', filepath);
end

[~, fname, fext] = fileparts(filepath);
fext = lower(fext);

switch fext
    case '.mat'
        datasets = readMatFile(filepath, fname);
    case '.csv'
        datasets = readCsvFile(filepath, fname, verbose);
    case '.xml'
        datasets = readXmlFile(filepath, fname, coordoffset);
    otherwise
        error('readCellLocationsFile:badFormat', ...
            ['Unsupported point file extension "%s" (%s). Supported: ' ...
             '.mat, .csv, .xml.'], fext, filepath);
end

%--------------------------------------------------------------------------
% resolve the channel of every dataset
%--------------------------------------------------------------------------
nds  = numel(datasets);
chan = chanopt;
if ~isempty(chan)
    chan = double(chan(:)).';
    if isscalar(chan)
        chan = repmat(chan, 1, nds);
    elseif numel(chan) ~= nds
        error('readCellLocationsFile:channelCount', ...
            ['''channel'' has %d entries but %s yielded %d point set(s). Pass ' ...
             'one channel per point set, or a single channel for all of them.'], ...
            numel(chan), filepath, nds);
    end
end

namechan = channelFromName(fname);
for k = 1:nds
    if ~isempty(chan)
        datasets(k).channel = chan(k);
        src = 'the ''channel'' option';
    elseif ~isnan(namechan)
        datasets(k).channel = namechan;
        src = 'the file name';
    elseif ~isnan(datasets(k).markertype)
        datasets(k).channel = datasets(k).markertype;
        src = 'the <Type> of the marker block';
    elseif ~isempty(defaultchan)
        datasets(k).channel = defaultchan;
        src = 'a fallback (NOT from the data)';
        warning('readCellLocationsFile:guessedChannel', ...
            ['Nothing identifies the channel of %s - its name carries no ' ...
             'chan<N> tag and the format holds none - so it is being counted ' ...
             'as channel %d. Pass ''channel'' to say which channel it really ' ...
             'is, or rename it to chan_<N>_%s.'], ...
            datasets(k).label, defaultchan, datasets(k).outbase);
    else
        error('readCellLocationsFile:unknownChannel', ...
            ['Cannot tell which channel %s belongs to. The file name carries ' ...
             'no chan<N> tag and the format holds no channel information. ' ...
             'Pass it explicitly, e.g. transformPointsToAtlas(..., ''channel'', 2).'], ...
            datasets(k).label);
    end

    % make sure the output files name their channel, without doubling up a
    % tag the source file already carried
    if isnan(channelFromName(datasets(k).outbase))
        datasets(k).outbase = sprintf('chan_%d_%s', datasets(k).channel, ...
            datasets(k).outbase);
    end

    if verbose
        fprintf('  %s -> channel %d (%d points, from %s)\n', ...
            datasets(k).label, datasets(k).channel, size(datasets(k).points, 1), src);
    end
end

end

%==========================================================================
% Format readers
%==========================================================================
function ds = readMatFile(filepath, fname)

dat = load(filepath);
if ~isfield(dat, 'cell_locations')
    error('readCellLocationsFile:noCellLocations', ...
        '%s does not contain a ''cell_locations'' variable.', filepath);
end

ds = emptyDataset();
ds.points     = single(dat.cell_locations);
ds.outbase    = regexprep(fname, '_sample$', '');
ds.sourcefile = filepath;
ds.sourcetype = 'mat';
ds.label      = [fname '.mat'];

if isfield(dat, 'cell_images') && ~isempty(dat.cell_images)
    ds.cell_images = dat.cell_images;
    ds.imwindow    = getOr(dat, 'imwindow', 6*[3 3 2]);
    if size(ds.cell_images, 1) ~= size(ds.points, 1)
        error('readCellLocationsFile:imageCountMismatch', ...
            ['%s holds %d cell_locations but %d cell_images rows; the file is ' ...
             'inconsistent and cannot be classified.'], filepath, ...
            size(ds.points, 1), size(ds.cell_images, 1));
    end
end

checkColumns(ds.points, filepath);

end

%--------------------------------------------------------------------------
function ds = readCsvFile(filepath, fname, verbose)

M = readmatrix(filepath);                 % skips a text header if present
M = M(~all(isnan(M), 2), :);              % drop blank lines

if isempty(M)
    error('readCellLocationsFile:emptyCsv', '%s holds no numeric rows.', filepath);
end

ds = emptyDataset();
ds.points     = single(M);
ds.outbase    = regexprep(fname, '_sample$', '');
ds.sourcefile = filepath;
ds.sourcetype = 'csv';
ds.label      = [fname '.csv'];

checkColumns(ds.points, filepath);

if verbose
    fprintf('  read %s: %d rows x %d columns, columns 1-3 taken as [x y z]\n', ...
        ds.label, size(M, 1), size(M, 2));
end

end

%--------------------------------------------------------------------------
function ds = readXmlFile(filepath, fname, coordoffset)

markers = readCellCounterXML(filepath, coordoffset);

ds = repmat(emptyDataset(), 1, numel(markers));
for k = 1:numel(markers)
    ds(k).points     = markers(k).points;
    ds(k).outbase    = sprintf('%s_type%d_cell_locations', fname, markers(k).type);
    ds(k).sourcefile = filepath;
    ds(k).sourcetype = 'xml';
    ds(k).markertype = markers(k).type;
    ds(k).label      = sprintf('%s.xml marker Type %d', fname, markers(k).type);
end

end

%==========================================================================
% Helpers
%==========================================================================
function ds = emptyDataset()
ds = struct('points', [], 'channel', NaN, 'outbase', '', 'sourcefile', '', ...
    'sourcetype', '', 'label', '', 'cell_images', [], 'imwindow', [], ...
    'markertype', NaN);
end

%--------------------------------------------------------------------------
function checkColumns(pts, filepath)
if size(pts, 2) < 3
    error('readCellLocationsFile:tooFewColumns', ...
        ['%s has %d column(s); at least 3 are needed ([x y z] in sample ' ...
         'space, optionally followed by descriptors).'], filepath, size(pts, 2));
end
end

%--------------------------------------------------------------------------
function ichan = channelFromName(fname)
%CHANNELFROMNAME Channel number encoded in a file name, or NaN.
%   Accepts the two conventions LightSuite writes ('chan_1_...' from the
%   lightsheet pipeline, 'chan03_...' from the slice module) plus a spelled-out
%   'channel3_' variant.
patterns = {'chan_(\d+)_', 'chan(\d+)_', 'channel_?(\d+)'};
for k = 1:numel(patterns)
    tok = regexp(fname, patterns{k}, 'tokens', 'once');
    if ~isempty(tok)
        ichan = str2double(tok{1});
        return;
    end
end
ichan = NaN;
end

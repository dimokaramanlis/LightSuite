function tfiles = resolveRawDataFiles(opts, datafolderoverride)
%RESOLVERAWDATAFILES Raw light-sheet TIFF list of a processed brain.
%   TFILES = RESOLVERAWDATAFILES(OPTS, DATAFOLDER) returns opts.tfiles when it
%   still points at existing files. Otherwise (raw data moved, or DATAFOLDER
%   given) it rescans DATAFOLDER, or opts.datafolder, exactly like
%   readLightsheetOpts does (natural sort of *.tif / *.tiff).
if nargin < 2; datafolderoverride = ''; end

tfiles = {};
if isempty(datafolderoverride) && isfield(opts, 'tfiles') && ~isempty(opts.tfiles) ...
        && exist(opts.tfiles{1}, 'file')
    tfiles = opts.tfiles;
    return;
end

datafolder = char(datafolderoverride);
if isempty(datafolder) && isfield(opts, 'datafolder')
    datafolder = opts.datafolder;
end
if isempty(datafolder) || ~exist(datafolder, 'dir')
    error('resolveRawDataFiles:notFound', ...
        ['Raw data not found. Stored tfiles are missing and datafolder ''%s'' ' ...
         'does not exist. Pass ''DataFolder'' to point at the relocated raw TIFFs.'], ...
        datafolder);
end
tiflist = cat(1, dir(fullfile(datafolder, '*.tif')), dir(fullfile(datafolder, '*.tiff')));
if isempty(tiflist)
    error('resolveRawDataFiles:notFound', 'No TIFF files found in %s.', datafolder);
end
tfiles = sort_nat(fullfile({tiflist(:).folder}', {tiflist(:).name}'));
end

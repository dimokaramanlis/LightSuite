function [volperm, trstruct, regopts, volname] = loadTracingRegistration(savepath, channel, volumepath)
% LOADTRACINGREGISTRATION  Load the volume and transform the tracing GUIs
%   (annotateGRINLens / annotateNeuropixelsProbes) annotate on.
%
%   [VOLPERM, TRSTRUCT, REGOPTS, VOLNAME] = loadTracingRegistration(SAVEPATH, CHANNEL, VOLUMEPATH)
%
%   Light-sheet brains: the 20 µm registration volume of the chosen channel
%   (resolveTracingVolumePath), permuted to atlas orientation with
%   permuteBrainVolume. Clicked points go through the full sample→atlas
%   transform (registrationPointsToAtlas).
%
%   Slice-module brains (isSliceRegistration): volume_registered/chanNN_*.tiff
%   is already the brain in Allen CCF space (10 µm, AP cropped to
%   atlasaplims, pages = ML). It is block-averaged to the registration
%   resolution (REGOPTS.registres) for display, giving [AP DV ML] voxels;
%   registrationPointsToAtlas then only undoes that downsampling and the AP
%   crop.
%
%   CHANNEL / VOLUMEPATH select the channel exactly as in
%   resolveTracingVolumePath; for slice brains CHANNEL indexes the
%   chanNN_*.tiff files and the default is the registration channel.
%
%   REGOPTS always has .registres. VOLNAME is the file name shown in the GUI.

    if nargin < 2; channel    = []; end
    if nargin < 3; volumepath = ''; end

    optsfile = fullfile(savepath, 'regopts.mat');
    trfile   = fullfile(savepath, 'transform_params.mat');
    if ~exist(optsfile, 'file') || ~exist(trfile, 'file')
        error('loadTracingRegistration:notFound', ...
            'regopts.mat or transform_params.mat not found in %s', savepath);
    end
    trstruct  = load(trfile);
    opts_data = load(optsfile);

    if ~isSliceRegistration(trstruct)
        regopts = opts_data.opts;
        volpath = resolveTracingVolumePath(savepath, regopts, channel, volumepath);
        fprintf('Loading volume from %s ...\n', volpath);
        volperm = permuteBrainVolume(readDownStack(volpath), trstruct.how_to_perm);
    else
        regopts = opts_data;                    % slice regopts are saved unwrapped
        volpath = resolveSliceVolumePath(savepath, channel, volumepath);
        ds      = regopts.registres / trstruct.atlasres;
        fprintf('Loading slice-registered volume from %s ...\n', volpath);
        volperm = readDownsampledStack(volpath, ds);
    end
    [~, n, e] = fileparts(volpath);
    volname   = [n e];
    fprintf('Volume size (permuted): %d x %d x %d\n', size(volperm));
end

%==========================================================================

function volpath = resolveSliceVolumePath(savepath, channel, volumepath)
    voldir = fullfile(savepath, 'volume_registered');
    if ~isempty(volumepath)
        volumepath = char(volumepath);
        volpath    = volumepath;
        if ~exist(volpath, 'file')
            [~, n, e] = fileparts(volumepath);
            volpath   = fullfile(voldir, [n e]);
        end
    else
        flist = [dir(fullfile(voldir, 'chan*.tif')); dir(fullfile(voldir, 'chan*.tiff'))];
        if isempty(flist)
            error('loadTracingRegistration:notFound', ...
                ['No chanNN_*.tif(f) in %s.\n' ...
                 'Run generateRegisteredSliceVolume first.'], voldir);
        end
        nums  = cellfun(@(s) str2double(regexp(s, '^chan(\d+)', 'tokens', 'once')), {flist.name});
        names = cellfun(@(s) regexprep(s, '^chan\d+_|\.tiff?$', ''), {flist.name}, 'UniformOutput', false);
        if isempty(channel)
            channel = defaultSliceChannel(savepath, names, nums);
        end
        k = find(nums == channel, 1);
        if isempty(k)
            error('loadTracingRegistration:notFound', ...
                'Channel %d not found in %s (available: %s).', channel, voldir, num2str(nums));
        end
        volpath = fullfile(voldir, flist(k).name);
    end
    if ~exist(volpath, 'file')
        error('loadTracingRegistration:notFound', 'Tracing volume not found:\n  %s', volpath);
    end
end

function channel = defaultSliceChannel(savepath, names, nums)
% registration channel from sliceinfo.regchan, else the first channel
    channel = min(nums);
    sfile   = fullfile(savepath, 'sliceinfo.mat');
    if exist(sfile, 'file')
        s = load(sfile, 'sliceinfo');
        regchan = getOr(getOr(s, 'sliceinfo', struct), 'regchan', '');
        k = find(strcmpi(names, regchan), 1);
        if ~isempty(k); channel = nums(k); end
    end
end

function vol = readDownsampledStack(volpath, ds)
% block-mean a multipage TIFF by integer factor DS in all three dims, reading
% one page at a time; trailing voxels that do not fill a block are dropped
    info = imfinfo(volpath);
    sz   = floor([info(1).Height, info(1).Width, numel(info)] / ds);
    vol  = zeros(sz, 'single');
    rows = 1:sz(1)*ds;  cols = 1:sz(2)*ds;
    for k = 1:sz(3)
        acc = zeros(sz(1)*ds, sz(2)*ds, 'single');
        for j = (k-1)*ds + (1:ds)
            im  = imread(volpath, 'Index', j, 'Info', info);
            acc = acc + single(im(rows, cols));
        end
        vol(:, :, k) = squeeze(mean(mean(reshape(acc, ds, sz(1), ds, sz(2)), 1), 3)) / ds;
    end
    vol = uint16(vol);
end

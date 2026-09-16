function [atlasvol, areastats] = generateRegisteredCordVolume(savepath, varargin)
%GENERATEREGISTEREDCORDVOLUME Put every channel of a cord sample into atlas space.
%
%   ATLASVOL = GENERATEREGISTEREDCORDVOLUME(SAVEPATH) replays the registration
%   on the image data itself and returns the result *in atlas space*: an array
%   of size [atlas volume, Nchannels] on the native spinal cord atlas grid, so
%   the annotation from LOADSPINALCORDATLAS indexes it directly, voxel for
%   voxel. It is the cord counterpart of GENERATEREGISTEREDBRAINVOLUMES.
%
%   Each channel goes through the same chain as a point does in CORDPOINTSTOATLAS
%   - permute, crop, straighten slice by slice, B-spline, affine, undo the
%   rostrocaudal flip - and is then resampled onto the native atlas grid. The
%   volumes are read back from the downsampled TIFFs written by
%   PREPROCESSLIGHTSHEETVOLUME, not from the raw data.
%
%   [ATLASVOL, AREASTATS] = GENERATEREGISTEREDCORDVOLUME(...) also returns the
%   per-region, per-segment summary of every channel, Nareas x Nsegments x
%   Nchannels. It is saved per channel as 'chanNN_intensities.mat' inside
%   SAVEPATH/volume_registered, next to the registered volumes.
%
%   Optional name-value arguments
%     'saveregisteredvolume' - write the registered volumes to
%                              SAVEPATH/volume_registered (default true).
%     'computestats'         - compute the per-region statistics (default true).
%     'writetocsv'           - also write the statistics as CSV (default from
%                              opts.writetocsv, else false).
%     'areafun'              - per-region summary statistic (default @median);
%                              any handle taking a vector and returning a scalar.
%     'outputres'            - 'atlas' (default) returns the native atlas grid;
%                              'registration' stops at the isotropic
%                              registration grid, which is smaller and is what
%                              the registration itself was computed on.
%
%   See also GENERATEREGISTEREDBRAINVOLUMES, TRANSFORMCORDPOINTSTOATLAS,
%   MULTIOBJCORDREGISTRATION, LOADSPINALCORDATLAS.

%==========================================================================
% 1. configuration
%==========================================================================
opts     = loadRegOpts(savepath);
savepath = opts.savepath;
trstruct = load(fullfile(savepath, 'transform_params.mat'));

defaultWriteCsv = getOr(opts, 'writetocsv', false);

p = inputParser;
addParameter(p, 'saveregisteredvolume', true, @(x) islogical(x) || isscalar(x));
addParameter(p, 'computestats',         true, @(x) islogical(x) || isscalar(x));
addParameter(p, 'writetocsv', defaultWriteCsv, @(x) islogical(x) || isscalar(x));
addParameter(p, 'areafun',          @median, @(x) isa(x, 'function_handle'));
addParameter(p, 'outputres',        'atlas', @(x) ischar(x) || isstring(x));
parse(p, varargin{:});
params    = p.Results;
outputres = validatestring(params.outputres, {'atlas', 'registration'});
%==========================================================================
registerpath = fullfile(savepath, 'volume_registered');
makeNewDir(registerpath);

Nchannels = opts.Nchans;
channames = getOr(opts, 'channames', repmat({''}, Nchannels, 1));
%==========================================================================
% 2. replay the registration on every channel
%==========================================================================
fprintf('Applying transforms... \n'); savetic = tic;

tforms     = trstruct.slicetforms;
sizetv     = trstruct.atlassize;
raout      = imref2d(sizetv([1 2]));
Rmoving    = imref3d([sizetv([1 2]) numel(tforms)]);
Rfixed     = imref3d(sizetv);

yrange = trstruct.samp_ikeepy;
xrange = trstruct.samp_ikeepx;
zrange = trstruct.samp_ikeeplong;

regvol = zeros([sizetv Nchannels], 'uint16');
for ichan = 1:Nchannels
    volpath = dir(fullfile(savepath, sprintf('chan_%d_*register*.tif', ichan)));
    assert(~isempty(volpath), 'generateRegisteredCordVolume:noChannelVolume', ...
        'No downsampled volume found for channel %d in %s.', ichan, savepath);
    currfname = fullfile(volpath(1).folder, volpath(1).name);
    fprintf('Registering %s\n', currfname);

    currvol = readDownStack(currfname);
    currvol = permute(currvol, trstruct.how_to_perm);
    currvol = currvol(yrange(1):yrange(2), xrange(1):xrange(2), zrange(1):zrange(2));

    assert(size(currvol, 3) == numel(tforms), ...
        'generateRegisteredCordVolume:sliceMismatch', ...
        ['Channel %d has %d slices in the kept range but the straightening ' ...
         'covers %d. The registration and the volumes are out of sync; re-run ' ...
         'the initialisation.'], ichan, size(currvol, 3), numel(tforms));

    currvol   = transformCordImageSlices(currvol, tforms, raout);
    volumereg = transformix(currvol, trstruct.tform_bspline_samp20um_to_atlas_20um_px, ...
        'movingscale', trstruct.registrationres(1) * 1e-3 * [1 1 1]);
    volumereg = uint16(abs(volumereg));

    regvol(:, :, :, ichan) = imwarp(volumereg, Rmoving, ...
        trstruct.tform_affine_samp20um_to_atlas_20um_px, 'OutputView', Rfixed);
    fprintf('Channel %d/%d done. Time %2.2f s. \n', ichan, Nchannels, toc(savetic));
end
%==========================================================================
% the atlas was flipped to match the sample during registration; put the
% result back into the atlas as it is distributed
if trstruct.tofliprc
    regvol = flip(regvol, 3);
end
%==========================================================================
% 3. onto the requested grid
%==========================================================================
[av, parcelinfo, segmentinfo, atlasres, nativesize] = loadCordAnnotationForOutput(trstruct, outputres);

if strcmp(outputres, 'atlas') && ~isequal(size(regvol, 1:3), nativesize)
    fprintf('Resampling to the native atlas grid (%s)... ', mat2str(nativesize)); tic;
    atlasvol = zeros([nativesize Nchannels], 'uint16');
    for ichan = 1:Nchannels
        atlasvol(:, :, :, ichan) = imresize3(regvol(:, :, :, ichan), nativesize);
    end
    fprintf('Done! Took %2.2f s.\n', toc);
else
    atlasvol = regvol;
end
clear regvol;
%==========================================================================
% 4. save
%==========================================================================
if params.saveregisteredvolume
    fprintf('Saving registered volumes... '); savetic = tic;
    saveLargeSliceVolume(permute(atlasvol, [1 2 4 3]), channames, registerpath);
    fprintf('Done! Took %2.2f s. \n', toc(savetic));
end
%==========================================================================
% 5. per-region, per-segment statistics
%==========================================================================
areastats = [];
if ~params.computestats
    return
end

fprintf('Summarising each region and segment with %s...\n', func2str(params.areafun));
proctic  = tic;
grouping = cordAtlasGrouping(av, segmentinfo, nativesize(3));

% names of every region at all three levels, so the outputs identify a region
% the way the brain ones do rather than by atlas id alone
areahierarchy = cordAreaHierarchy(parcelinfo, grouping.avinds);

areastats = nan(grouping.Nareas, grouping.Nsegments, Nchannels, 'single');
for ichan = 1:Nchannels
    [medianoverareas, volumeoverareas] = cordAreaStatistics(...
        atlasvol(:, :, :, ichan), av, grouping, params.areafun, atlasres);

    areastats(:, :, ichan) = medianoverareas;

    areaidx     = grouping.avinds;
    segmentname = grouping.segnames;
    areafun     = params.areafun;
    areafunname = func2str(areafun);
    fmatname    = fullfile(registerpath, sprintf('chan%02d_intensities.mat', ichan));
    save(fmatname, 'medianoverareas', 'areaidx', 'volumeoverareas', ...
        'segmentname', 'areahierarchy', 'areafun', 'areafunname');

    if params.writetocsv
        currtable = cordAreaTable(areahierarchy, segmentname, ...
            struct('Intensity', medianoverareas, 'Volume_mm3', volumeoverareas));
        writetable(currtable, ...
            fullfile(registerpath, sprintf('chan%02d_intensities.csv', ichan)), 'Delimiter', ';');
    end
    fprintf('Channel %d/%d done. Time %2.2f s. \n', ichan, Nchannels, toc(proctic));
end
%==========================================================================
end

%==========================================================================
% Local helpers
%==========================================================================
function [av, parcelinfo, segmentinfo, atlasres, nativesize] = loadCordAnnotationForOutput(trstruct, outputres)
%LOADCORDANNOTATIONFOROUTPUT Annotation on the grid the output volume lives on.
%   Never flipped: the registered volume has already been flipped back, so the
%   atlas is needed exactly as distributed.

switch outputres
    case 'atlas'
        [~, av, parcelinfo, segmentinfo, atlasres, nativesize] = loadSpinalCordAtlas();
    case 'registration'
        registrationres = trstruct.registrationres(1) * [1 1 1];
        [~, av, parcelinfo, segmentinfo, ~, nativesize] = loadSpinalCordAtlas(registrationres);
        atlasres = registrationres;
end

end

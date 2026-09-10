function opts = prepareCordSampleForRegistration(inputpath, varargin)
%PREPARECORDSAMPLEFORREGISTRATION First half of the spinal cord initialisation.
%
%   OPTS = PREPARECORDSAMPLEFORREGISTRATION(INPUTPATH) takes the downsampled
%   registration volume written by PREPROCESSLIGHTSHEETVOLUME and turns it into
%   something the aligner can work on. INPUTPATH is the folder holding
%   regopts.mat (OPTS.savepath); the options struct is updated in place and
%   saved back.
%
%   This is where a cord sample departs from a brain. A brain is a compact
%   object that only has to be re-oriented, so INITIALIZEREGISTRATION can go
%   straight to a similarity transform. A cord is a long, bent, twisted tube cut
%   out of a larger block, so before any atlas can be fitted the pipeline has to
%   work out how it is lying in the volume:
%
%     1. Which array dimension is the long (rostrocaudal) axis, and move it last.
%     2. Segment the cord: Otsu threshold, dilate, keep the largest connected
%        component. The cross-sectional area of that component per slice drives
%        the next two steps.
%     3. Which end is rostral. If the sample runs caudorostral the *atlas* is
%        flipped to match, not the sample - see LOADCORDATLASVOLUMES.
%     4. Where the cord stops and brain tissue starts. Slices whose area is far
%        above the rest of the cord are dropped, together with everything beyond
%        them; the kept range is OPTS.ikeeprange.
%     5. Crop the remaining bounding box so the volume is not mostly background.
%     6. Predict, for every remaining slice, the centre of the section and its
%        dorsoventral angle (AUTOCORDCENTERLINE). SPINAL_CORD_ALIGNER fits these
%        predictions together with whatever the user clicks, dropping them
%        wherever a user point is close enough to speak for that stretch of cord
%        ('userexclusionradius' below).
%
%   The cropped registration volume is written next to the other volumes as
%   'cord_sample_register_<res>um.tif' and referenced by OPTS.cordvolpath, so
%   nothing large ends up inside regopts.mat.
%
%   Run SPINAL_CORD_ALIGNER next, then INITIALIZECORDREGISTRATION.
%
%   Optional name-value arguments
%     'Volume'           - use this volume instead of reading OPTS.regvolpath.
%     'longaxis'         - force the long axis (1, 2 or 3) instead of taking the
%                          largest dimension.
%     'cropmargin'       - pixels of background kept around the cord (default 2).
%     'brainthreshold'   - how many robust standard deviations above the median
%                          cross-section still counts as cord (default 3).
%     'centerlinelambda' - smoothing of the automatic centre line (default 25).
%     'userexclusionradius' - how far a click in SPINAL_CORD_ALIGNER reaches, in
%                          slices (default 100). Inside that distance of a user
%                          point the automatic prediction is dropped from the
%                          fit, so the user's points decide that stretch of cord
%                          on their own. Raise it for a coarsely sampled cord or
%                          when corrections keep being pulled back towards the
%                          prediction; it can also be changed inside the GUI.
%     'extractclouds'    - also extract sample and atlas point clouds, for
%                          point-cloud based initialisation (default false; the
%                          image-based path used by the pipeline does not need
%                          them).
%
%   See also PREPROCESSLIGHTSHEETVOLUME, SPINAL_CORD_ALIGNER,
%   INITIALIZECORDREGISTRATION, AUTOCORDCENTERLINE.

%==========================================================================
p = inputParser;
addRequired(p,  'inputpath', @(x) ischar(x) || isstring(x) || isstruct(x));
addParameter(p, 'Volume',           [], @isnumeric);
addParameter(p, 'longaxis',         [], @(x) isempty(x) || (isscalar(x) && ismember(x, 1:3)));
addParameter(p, 'cropmargin',        2, @(x) isscalar(x) && x >= 0);
addParameter(p, 'brainthreshold',    3, @(x) isscalar(x) && x > 0);
addParameter(p, 'centerlinelambda', 25, @(x) isscalar(x) && x >= 0);
addParameter(p, 'userexclusionradius', 100, @(x) isscalar(x) && x >= 0);
addParameter(p, 'extractclouds', false, @(x) islogical(x) || isscalar(x));
parse(p, inputpath, varargin{:});
params = p.Results;
%==========================================================================
opts                     = loadRegOpts(inputpath);
opts.samplekind          = 'cord';
opts.registrationres     = opts.registres * [1 1 1];
opts.userexclusionradius = round(params.userexclusionradius);
%==========================================================================
if ~isempty(params.Volume)
    regvol = params.Volume;
else
    assert(isfield(opts, 'regvolpath') && exist(opts.regvolpath, 'file') == 2, ...
        'prepareCordSampleForRegistration:noVolume', ...
        ['No registration volume found. Run preprocessLightSheetVolume first, ' ...
         'or pass one with the ''Volume'' option.']);
    fprintf('Reading registration volume %s\n', opts.regvolpath);
    regvol = readDownStack(opts.regvolpath);
end
opts.regvolsize = size(regvol, 1:3);
%==========================================================================
% 1. the long axis goes last
if isempty(params.longaxis)
    [~, ilong] = max(opts.regvolsize);
else
    ilong = params.longaxis;
end
axperm          = setdiff(1:3, ilong);
opts.sampleperm = [axperm ilong];
regvol          = permute(regvol, opts.sampleperm);
regvolsize      = size(regvol, 1:3);
fprintf('Long axis is dimension %d of the sample, moved to last (%s).\n', ...
    ilong, mat2str(opts.sampleperm));
%==========================================================================
% 2. segment the cord
fprintf('Segmenting the cord... '); tic;
isamprand           = randperm(numel(regvol), min(2e4, numel(regvol)));
sampsuse            = regvol(isamprand);
backval             = mode(sampsuse(sampsuse > 0));
regvol(regvol == 0) = backval;

binvol    = imbinarize(regvol);
binvol    = imdilate(binvol, strel('cuboid', [3 3 3]));
cc        = bwconncomp(binvol, 18);
cinfo     = regionprops3(cc, 'Volume', 'VoxelList');
[~, imax] = max(cinfo.Volume);
voxset    = cinfo.VoxelList{imax};                   % columns are [x y z]
cordmask  = false(regvolsize);
cordmask(sub2ind(regvolsize, voxset(:, 2), voxset(:, 1), voxset(:, 3))) = true;
fprintf('Done! Took %2.2f s.\n', toc);
%==========================================================================
% 3. rostrocaudal direction, from how the cord thickens towards the brain
cordarea      = squeeze(sum(cordmask, 1:2));
Ncheck        = ceil(0.05 * numel(cordarea));
opts.tofliprc = mean(cordarea(1:Ncheck)) < mean(cordarea(end-Ncheck+1:end));
if opts.tofliprc
    fprintf('Sample runs caudorostral; the atlas will be flipped to match it.\n');
end
opts.cordarea = cordarea;
%==========================================================================
% 4. drop the brain end, if part of it was imaged
ihigh = cordarea > (median(cordarea) + params.brainthreshold * robustStd(cordarea));
fprintf('%2.2f%% of the slices are wider than a cord should be.\n', mean(ihigh) * 100);
ikeep = [1 regvolsize(3)];
if mean(ihigh) > 0.01
    if opts.tofliprc
        ilast = find(ihigh == 1, 1);
        ikeep = [1 ilast + 1];
    else
        ifirst = find(ihigh == 0, 1);
        ikeep  = [ifirst - 1 regvolsize(3)];
    end
    ikeep = [max(ikeep(1), 1) min(ikeep(2), regvolsize(3))];
    fprintf('Keeping slices %d-%d of %d along the long axis.\n', ...
        ikeep(1), ikeep(2), regvolsize(3));
end
opts.ikeeprange = ikeep;
%==========================================================================
% 5. crop to the cord, with a small margin
mrg         = params.cropmargin;
xrange      = [max(min(voxset(:, 1)) - mrg, 1) min(max(voxset(:, 1)) + mrg, regvolsize(2))];
yrange      = [max(min(voxset(:, 2)) - mrg, 1) min(max(voxset(:, 2)) + mrg, regvolsize(1))];
opts.xrange = xrange;
opts.yrange = yrange;

cropvol  = regvol(  yrange(1):yrange(2), xrange(1):xrange(2), ikeep(1):ikeep(2));
cropmask = cordmask(yrange(1):yrange(2), xrange(1):xrange(2), ikeep(1):ikeep(2));
fprintf('Cropped the registration volume to %s px.\n', mat2str(size(cropvol)));
%==========================================================================
% 6. predict the centre line and the twist
fprintf('Predicting the cord centre line...\n'); tic;
opts.cordauto = autoCordCenterline(cropvol, cropmask, 'lambda', params.centerlinelambda);
fprintf('Done! Took %2.2f s.\n', toc);

cf = plotCordCenterline(cropvol, cropmask, opts.cordauto);
print(cf, fullfile(opts.savepath, 'cord_automatic_centerline'), '-dpng');
close(cf);
%==========================================================================
% save the cropped volume the aligner and the straightening work on
opts.cordvolpath = fullfile(opts.savepath, ...
    sprintf('cord_sample_register_%dum.tif', opts.registres));
saveopts = struct('compress', 'lzw', 'message', false);
if exist(opts.cordvolpath, 'file')
    delete(opts.cordvolpath);
end
saveastiff(cropvol, opts.cordvolpath, saveopts);
%==========================================================================
% optional point clouds, for point-cloud based initialisation
if params.extractclouds
    fprintf('Extracting point clouds... '); tic;
    opts.smpts    = spinalCordPointCloud(cropvol, [1 size(cropvol, 3)], cropmask);
    [~, ~, tvpts] = loadSpinalCordAtlasAndPoints(opts.registrationres);
    if opts.tofliprc
        [~, ~, atlasinfo] = loadCordAtlasVolumes(opts);
        tvpts(:, 3) = atlasinfo.regsize(3) + 1 - tvpts(:, 3);
    end
    opts.tvpts = tvpts;
    fprintf('Done! Took %2.2f s.\n', toc);
end
%==========================================================================
saveRegOpts(opts);
fprintf(['Sample prepared. Refine the centre line with spinal_cord_aligner, ' ...
    'then run initializeCordRegistration.\n']);
%==========================================================================
end

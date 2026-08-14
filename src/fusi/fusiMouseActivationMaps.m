function res = fusiMouseActivationMaps(pathrawdata, mousename, opts)
%FUSIMOUSEACTIVATIONMAPS Per-mouse fUS activation maps and timecourses.
%   res = fusiMouseActivationMaps(pathrawdata, mousename) runs the whole
%   activation pipeline for one mouse: it computes per-session maps and
%   peristimulus timecourses (fusiSessionActivationMaps), then aggregates
%   them across sessions. When the mouse has been registered to the Allen
%   atlas (regopts.mat + transform_params.mat present), maps and per-area,
%   block-triggered timecourses are moved to atlas space with
%   warpFusiSessionToAtlas / applyFusiTransforms and averaged there; otherwise
%   maps are averaged in the mouse's own seed space via the anatomy transforms.
%   In atlas space the results are kept for BOTH the full non-rigid warp
%   (atlasMaps / areaMapVals / areaTC) and the affine-only warp (atlasMapsAff /
%   areaMapValsAff / areaTCaff), so the two registrations can be compared.
%
%   Input:
%     pathrawdata  root '...\fus\anes'.
%     mousename    e.g. 'DS_WT61'.
%     opts (optional):
%        .finshape        native volume shape (default [36 64 54]).
%        .space           'auto' (default), 'atlas' or 'seed'.
%        .areatimecourses per-area atlas timecourses too (default true).
%        .statclim        montage colour limits (default [-8 8]).
%        .plot            save figures (default true).
%        .maxsessions     cap the number of sessions (default inf; for quick
%                         runs / testing).
%        .verbose         default true.
%        .sessopts        struct forwarded to fusiSessionActivationMaps.
%
%   Output struct res (also saved to <lightsuite>/activation_maps/
%   <mouse>_mouse_actmaps.mat): see fields set below.
%
%   See also FUSISESSIONACTIVATIONMAPS, WARPFUSISESSIONTOATLAS,
%   FUSIPERISTIMULUSTIMECOURSES, PLOTFUSIACTIVATIONMONTAGE, PLOTFUSITIMECOURSES.

if nargin < 3, opts = struct(); end
opts = setdefault(opts, 'finshape', [36 64 54]);
opts = setdefault(opts, 'space', 'auto');
opts = setdefault(opts, 'areatimecourses', true);
opts = setdefault(opts, 'statclim', [-8 8]);
opts = setdefault(opts, 'plot', true);
opts = setdefault(opts, 'verbose', true);
opts = setdefault(opts, 'maxsessions', inf);   % cap sessions (quick runs)
opts = setdefault(opts, 'sessopts', struct());

savepath = fullfile(pathrawdata, mousename, 'lightsuite');
mapdir   = fullfile(savepath, 'activation_maps');
makeNewDir(mapdir);

% ---- discover sessions -------------------------------------------------
ap = dir(fullfile(pathrawdata, mousename, '**', '*_FUS.mat'));
ap = ap(~contains({ap.name}, 'template'));
sesspaths    = fullfile({ap.folder}', {ap.name}');
sessionnames = erase({ap.name}', '_FUS.mat');
Nsess        = numel(sesspaths);
if isfinite(opts.maxsessions) && opts.maxsessions < Nsess
    Nsess = opts.maxsessions;
    sesspaths = sesspaths(1:Nsess); sessionnames = sessionnames(1:Nsess);
end
if opts.verbose, fprintf('\n##### %s: %d sessions #####\n', mousename, Nsess); end

% ---- anatomy (session voxel size, + seed-space transforms) -------------
anatfile = fullfile(savepath, sprintf('%s_anatomy.mat', mousename));
anat = [];
if exist(anatfile, 'file')
    anat = load(anatfile, 'tforms', 'session_names', 'session_voxelsize_mm', ...
        'voxelsize_mm', 'volumeavg');   % skip the big sessionvolumes array
end
voxsize = [0.1971 0.15 0.15];
if ~isempty(anat) && isfield(anat, 'session_voxelsize_mm')
    voxsize = anat.session_voxelsize_mm;
end

% ---- decide output space -----------------------------------------------
haveAtlas = exist(fullfile(savepath, 'regopts.mat'), 'file') && ...
            exist(fullfile(savepath, 'transform_params.mat'), 'file');
useAtlas = (strcmp(opts.space, 'atlas')) || (strcmp(opts.space, 'auto') && haveAtlas);
if useAtlas && ~haveAtlas
    error('fusiMouseActivationMaps:noAtlas', ...
        'space=atlas requested but regopts/transform_params missing for %s', mousename);
end
optA = [];
if useAtlas
    optA = load(fullfile(savepath, 'regopts.mat'));
    if isempty(anat) || ~isfield(anat, 'tforms')
        error('fusiMouseActivationMaps:noTforms', ...
            'atlas warp needs %s_anatomy.mat with .tforms (session->seed).', mousename);
    end
end

% ---- per-session maps/timecourses, aggregated on the fly ---------------
sessopts = mergestruct(struct('finshape', opts.finshape, 'timecourses', true, ...
    'savepath', mapdir, 'verbose', opts.verbose), opts.sessopts);

roi = struct('combined', [], 'object', [], 'scrambled', []);   % [nlags x Nsess]
lags = []; mapnames = {}; groupidx = [];
if useAtlas
    asz = size(optA.atlas);
    mapSum    = zeros([asz 4], 'single'); mapCnt    = zeros([asz 4], 'single');
    mapSumAff = zeros([asz 4], 'single'); mapCntAff = zeros([asz 4], 'single');
    areaValsAcc = []; areaTCAcc = []; areaValsAccAff = []; areaTCAccAff = [];
else
    Rseed  = imref3d(size(anat.volumeavg));
    facvol = voxsize ./ anat.voxelsize_mm;
    seedSum = struct(); seedCnt = struct();
    for f = {'combined','object','scrambled'}
        seedSum.(f{1}) = zeros([size(anat.volumeavg) 1], 'single');
        seedCnt.(f{1}) = zeros([size(anat.volumeavg) 1], 'single');
    end
end

for i = 1:Nsess
    sr = fusiSessionActivationMaps(sesspaths{i}, sessopts);
    % ROI (native) peristimulus traces, averaged across sessions later
    if ~isempty(sr.timecourses)
        lags = sr.timecourses.lags;
        for f = {'combined','object','scrambled'}
            roi.(f{1})(:, i) = sr.timecourses.roi.(f{1}).mean;
        end
    end
    if useAtlas
        % session->seed rigid transform from the anatomy (no re-fitting)
        ia = find(strcmp(anat.session_names, sessionnames{i}), 1);
        if isempty(ia)
            warning('%s: session %s not in anatomy tforms; skipping atlas warp.', ...
                mousename, sessionnames{i});
            continue;
        end
        atl = warpFusiSessionToAtlas(optA, anat.tforms{ia}, voxsize, sr, ...
            struct('finshape', opts.finshape, ...
            'areatimecourses', opts.areatimecourses, 'verbose', opts.verbose));
        stack = @(m) cat(4, m.combined_corr, m.object_corr, m.scrambled_corr, m.combined_tscore);
        av  = stack(atl.maps);    valid  = isfinite(av);  av(~valid)   = 0;
        avA = stack(atl.mapsAff); validA = isfinite(avA); avA(~validA) = 0;
        mapSum    = mapSum    + av;  mapCnt    = mapCnt    + valid;
        mapSumAff = mapSumAff + avA; mapCntAff = mapCntAff + validA;
        mapnames = atl.mapnames; groupidx = atl.groupidx;
        if isempty(areaValsAcc)
            areaValsAcc    = nan([size(atl.mapAreaVals) Nsess], 'single');
            areaValsAccAff = nan([size(atl.mapAreaVals) Nsess], 'single');
        end
        areaValsAcc(:, :, i)    = atl.mapAreaVals;
        areaValsAccAff(:, :, i) = atl.mapAreaValsAff;
        if ~isempty(atl.areaTC)
            if isempty(areaTCAcc)
                areaTCAcc    = nan([size(atl.areaTC) Nsess], 'single');
                areaTCAccAff = nan([size(atl.areaTC) Nsess], 'single');
            end
            areaTCAcc(:, :, :, i)    = atl.areaTC;
            areaTCAccAff(:, :, :, i) = atl.areaTCaff;
        end
    else
        ia = find(strcmp(anat.session_names, sessionnames{i}), 1);
        if isempty(ia), continue; end
        for f = {'combined','object','scrambled'}
            w = warpToSeed(sr.maps.(f{1}).corr, facvol, anat.tforms{ia}, Rseed);
            v = isfinite(w); w(~v) = 0;
            seedSum.(f{1}) = seedSum.(f{1}) + w;
            seedCnt.(f{1}) = seedCnt.(f{1}) + v;
        end
    end
end

% ---- assemble results --------------------------------------------------
res = struct('mousename', mousename, 'Nsessions', Nsess, ...
    'sessionnames', {sessionnames}, 'lags', lags);
res.roiTC = struct();
for f = {'combined','object','scrambled'}
    R = roi.(f{1});
    res.roiTC.(f{1}) = struct('mean', mean(R, 2, 'omitnan'), ...
        'sem', std(R, 0, 2, 'omitnan') ./ sqrt(max(sum(~isnan(R), 2), 1)));
end

if useAtlas
    res.space = 'atlas';
    res.mapnames = mapnames; res.groupidx = groupidx;
    meanvol = @(S, C) setnan(S ./ max(C, 1), C == 0);
    meanmaps    = meanvol(mapSum,    mapCnt);
    meanmapsAff = meanvol(mapSumAff, mapCntAff);
    res.atlasMaps = struct(); res.atlasMapsAff = struct();
    for k = 1:numel(mapnames)
        res.atlasMaps.(mapnames{k})    = meanmaps(:, :, :, k);
        res.atlasMapsAff.(mapnames{k}) = meanmapsAff(:, :, :, k);
    end
    res.areaMapVals    = mean(areaValsAcc,    3, 'omitnan');  % [nArea x 4] non-rigid
    res.areaMapValsAff = mean(areaValsAccAff, 3, 'omitnan');  % [nArea x 4] affine
    if ~isempty(areaTCAcc)
        res.areaTC    = mean(areaTCAcc,    4, 'omitnan');     % [nArea x nlags x 3]
        res.areaTCaff = mean(areaTCAccAff, 4, 'omitnan');
    end
    res.atlas = optA.atlas;
else
    res.space = 'seed';
    res.seedMaps = struct();
    for f = {'combined','object','scrambled'}
        m = seedSum.(f{1}) ./ max(seedCnt.(f{1}), 1);
        m(seedCnt.(f{1}) == 0) = NaN;
        res.seedMaps.(f{1}) = m;
    end
    res.anatomy = anat.volumeavg;
end

save(fullfile(mapdir, sprintf('%s_mouse_actmaps.mat', mousename)), ...
    '-struct', 'res', '-v7.3');

% ---- figures -----------------------------------------------------------
if opts.plot
    plotMousefigures(res, mapdir, opts);
end
if opts.verbose, fprintf('##### %s done (%s space). #####\n', mousename, res.space); end
end

% =========================================================================
function plotMousefigures(res, mapdir, opts)
% combined-response montage over anatomy/atlas
if strcmp(res.space, 'atlas')
    under  = sqrt(single(res.atlas));
    sidx   = round(linspace(0.25, 0.80, 32) * size(res.atlas, 1));
    stat   = res.atlasMaps.combined_tscore;   % T-score
    saxis  = 1; clim = opts.statclim; thr = 3;
else
    under  = sqrt(max(res.anatomy, 0));
    sidx   = [];
    stat   = res.seedMaps.combined;           % correlation in seed space
    saxis  = 3; clim = [-0.5 0.5]; thr = 0;
end
mont = @(v, tag, ttl) plotFusiActivationMontage(v, struct('underlay', under, ...
    'clim', clim, 'thresh', thr, 'sliceaxis', saxis, 'sliceidx', sidx, ...
    'visible', 'off', 'title', ttl, 'savepng', ...
    fullfile(mapdir, sprintf('%s_mouse_%s', res.mousename, tag))));
mont(stat, 'combined', sprintf('%s combined visual response (%s non-rigid, N=%d)', ...
    res.mousename, res.space, res.Nsessions));
% affine-only version, for the non-rigid vs affine comparison
if strcmp(res.space, 'atlas')
    mont(res.atlasMapsAff.combined_tscore, 'combined_affine', ...
        sprintf('%s combined visual response (atlas affine-only, N=%d)', ...
        res.mousename, res.Nsessions));
end
% peristimulus ROI timecourses
if ~isempty(res.lags)
    plotFusiTimecourses(res.lags, res.roiTC, struct('visible', 'off', ...
        'title', sprintf('%s peristimulus response (responsive ROI, N=%d)', ...
        res.mousename, res.Nsessions), 'savepng', ...
        fullfile(mapdir, sprintf('%s_mouse_timecourses', res.mousename))));
end
end

% -------------------------------------------------------------------------
function w = warpToSeed(mapvol, facvol, tform, Rseed)
mapvol(~isfinite(mapvol)) = 0;
mapatl = imresize3(single(mapvol), 'Scale', facvol);
w = imwarp(mapatl, imref3d(size(mapatl)), tform, 'OutputView', Rseed, 'FillValues', nan);
end

function s = mergestruct(s, t)
if isempty(t) || ~isstruct(t), return; end
for f = fieldnames(t)', s.(f{1}) = t.(f{1}); end
end

function v = setnan(v, mask)
v(mask) = NaN;
end

function s = setdefault(s, f, v)
if ~isfield(s, f) || isempty(s.(f)), s.(f) = v; end
end

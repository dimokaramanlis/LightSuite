function res = buildFusiAnatomy(opts, sessionvols, voxelsize_mm, sessionnames, savesessions)
%BUILDFUSIANATOMY Anatomical template (atlas resolution) from repeated fUSI sessions.
%   res = buildFusiAnatomy(opts, sessionvols, voxelsize_mm) builds a within-mouse
%   anatomy scan out of many freehand-repositioned fUSI sessions by upsampling
%   every session to the atlas resolution and rigidly aligning them to a seed:
%     1) UPSAMPLE: each session is resampled once to isotropic atlas resolution
%        and equalized to a common median so no session dominates the average.
%     2) SEED: the user picks a seed session (selectSeedSession, cached in
%        <opts.savepath>/seed_session_for_anatomy.txt). The seed, on the atlas
%        grid, anchors the common "seed space".
%     3) REGISTRATION: every other session is rigidly registered to the seed
%        once with elastix, then warped into seed space. Each session-to-seed
%        transform is cached under <opts.savepath>/anatomy_registration/<session>/
%        and reloaded on reruns, so the elastix fit only runs the first time.
%     4) COMBINE: the anatomy is the mean over the aligned sessions (NaN-omitting)
%        and is written to <opts.savepath>/<opts.mousename>_anatomy.mat.
%
%   The template stays in the NATIVE session orientation (only the sampling
%   changes), so it feeds straight into setupFusiOptions, which decides the
%   sample-to-atlas permutation.
%
%   Inputs:
%     opts         - struct with fields savepath, mousename and atlas_res.
%     sessionvols  - [d1 d2 d3 Nsessions] per-session volumes (median over time
%                    of each recording). Coronal slices run along dim 3.
%     voxelsize_mm - [d1 d2 d3] voxel size of sessionvols, in mm.
%     sessionnames - (optional) cellstr of Nsessions names.
%     savesessions - (optional, default true) also store every session warped
%                    into seed space for QC.
%
%   Output (also the contents of <mousename>_anatomy.mat):
%     res.volume         - [D1 D2 D3] single anatomy at atlas resolution,
%                          ready to feed into setupFusiOptions.
%     res.voxelsize_mm   - resolution of res.volume (== opts.atlas_res)
%     res.seed_session   - index of the seed session
%     res.session_names  - names of all sessions, in order
%     res.session_voxelsize_mm - resolution of the input volumes
%     res.tforms         - Nsessions cell of session->seed rigid transforms
%     res.sessionvolumes - [D1 D2 D3 Nsessions] single sessions warped into seed
%                          space at atlas resolution (QC only, optional)

Nsessions = size(sessionvols, 4);

if nargin < 4 || isempty(sessionnames)
    sessionnames = arrayfun(@(x) sprintf('session%02d', x), 1:Nsessions, 'UniformOutput', false);
end
if nargin < 5 || isempty(savesessions), savesessions = true; end
voxelsize_mm = double(voxelsize_mm(:)).';
atlasres     = double(opts.atlas_res(:)).';
if isscalar(atlasres), atlasres = atlasres([1 1 1]); end

makeNewDir(opts.savepath);
%==========================================================================
% (1) upsample every session to isotropic atlas resolution once, cleaned and
% equalized to a common median so no session dominates the average
facvol    = voxelsize_mm./atlasres;
sessatlas = cell(Nsessions, 1);
for isess = 1:Nsessions
    v = sanitizeVol(imresize3(single(sessionvols(:, :, :, isess)), 'Scale', facvol));
    m = median(v(v > 0), 'omitnan');
    if ~isfinite(m) || m <= 0, m = 1; end
    sessatlas{isess} = v./m;
end
%==========================================================================
% (2) pick the seed session that defines the anatomy space
iseed   = selectSeedSession(sessionvols, opts.savepath, sessionnames, voxelsize_mm);
seedvol = sessatlas{iseed};
Rseed   = imref3d(size(seedvol));
fprintf('Building anatomy for %s: %d sessions, seed %d (%s), atlas res %.0f um, rigid transform.\n', ...
    opts.mousename, Nsessions, iseed, sessionnames{iseed}, atlasres(1)*1e3);
%==========================================================================
% (3) rigidly register each session to the seed once (cached) and warp it into
% seed space
regpath = fullfile(opts.savepath, 'anatomy_registration');
makeNewDir(regpath);
% coarse-to-fine elastix pyramid to align a freehand session from scratch
% elastixparams = struct('NumberOfResolutions', 2, 'ImagePyramidSchedule', [2 1], ...
%     'MaximumNumberOfIterations', [1000 1000]);

elastixparams = struct('NumberOfResolutions', 2, 'ImagePyramidSchedule', [2 1], ...
    'MaximumNumberOfIterations', [1000 1000], 'NumberOfHistogramBins', 64);


tforms  = cell(Nsessions, 1);
aligned = nan([size(seedvol) Nsessions], 'single');   % each session in seed space
for isess = 1:Nsessions
    outpath = fullfile(regpath, sessionnames{isess});
    makeNewDir(outpath);
    if isess == iseed
        tforms{isess}           = rigidtform3d(eye(4));   % seed anchors the frame
        aligned(:, :, :, isess) = seedvol;
        continue;
    end
    volatl  = sessatlas{isess};
    tform   = registerSessionToSeed(seedvol, volatl, atlasres(1), outpath, elastixparams);
    warpvol = imwarp(volatl, imref3d(size(volatl)), tform, 'OutputView', Rseed, 'FillValues', nan);
    warpvol(warpvol < 0)    = 0;
    tforms{isess}           = tform;
    aligned(:, :, :, isess) = warpvol;
end
%==========================================================================
% (4) combine into the anatomy (mean over sessions, NaN-omitting)
volumeavg = mean(aligned, 4, 'omitnan');

res = struct();
res.volumeavg            = volumeavg;
res.voxelsize_mm         = atlasres;
res.seed_session         = iseed;
res.session_names        = sessionnames(:);
res.session_voxelsize_mm = voxelsize_mm;
res.tforms               = tforms;
if savesessions
    res.sessionvolumes = aligned;
end
%==========================================================================
% QC: overlay each session, warped into seed space, on the anatomy
for isess = 1:Nsessions
    cf = visualizeTransformMatch(forDisplay(sqrt(seedvol)), forDisplay(sqrt(aligned(:, :, :, isess))), 3);
    print(cf, fullfile(regpath, sessionnames{isess}, 'rigid_registration_to_template'), '-dpng')
    close(cf);
end
%==========================================================================
savename = fullfile(opts.savepath, sprintf('%s_anatomy.mat', opts.mousename));
save(savename, '-struct', 'res', '-v7.3');
fprintf('Anatomy volume [%s] at %.0f um saved to %s\n', ...
    num2str(size(res.volumeavg)), atlasres(1)*1e3, savename);
%==========================================================================
% one overview figure of the anatomy against the seed alone
cf = visualizeTransformMatch(forDisplay(sqrt(seedvol)), forDisplay(sqrt(volumeavg)), 3);
print(cf, fullfile(opts.savepath, sprintf('%s_anatomy', opts.mousename)), '-dpng')
close(cf);
%==========================================================================
end

function tform = registerSessionToSeed(seedvol, sessvol, volscale, outpath, elastixparams)
%REGISTERSESSIONTOSEED Elastix rigid fit of a session to the seed.
%   Registers the seed onto the session and inverts, so the returned rigid
%   transform maps the session into seed space. Cached per session (saved to
%   rigid_to_seed.txt and reloaded if it already exists) so reruns are cheap.
targetfile = fullfile(outpath, 'rigid_to_seed.txt');
if exist(targetfile, 'file')
    fprintf('    Found rigid transform in %s, loading...\n', outpath)
else
    fprintf('    Fitting rigid transform...\n')
    [~, ~, tformpath, ~] = performElastixRigidRegistration(seedvol, sessvol, volscale, outpath, elastixparams);
    movefile(tformpath, targetfile);
end
tform = parse_elastix_tform_all(targetfile);
tform = rigidtform3d(tform.A);
tform = tform.invert();
end

function vol = sanitizeVol(vol)
%SANITIZEVOL Force a volume to finite, non-negative single precision.
%   Real I./Inorm data can hold Inf/NaN (where Inorm == 0) and small negatives;
%   either makes elastix diverge or crash, so clean them before registration.
vol = single(vol);
vol(~isfinite(vol)) = 0;
vol(vol < 0)        = 0;
end

function im = forDisplay(vol)
%FORDISPLAY uint8 copy scaled to its 1-99% range, for the imfuse overlays in
%   visualizeTransformMatch. quantile ignores the NaN FOV border, and NaNs map
%   to 0 (black) in the uint8 cast, so no extra masking is needed.
vol    = single(vol);
imlims = quantile(vol, [0.01 0.99], 'all');
im     = uint8(255*(vol - imlims(1))/range(imlims));
end

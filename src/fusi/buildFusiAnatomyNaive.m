function res = buildFusiAnatomyNaive(opts, sessionvols, voxelsize_mm, sessionnames, savesessions, Npass, combinemode, transformtype)
%BUILDFUSIANATOMYNAIVE Atlas-resolution anatomical template from repeated fUSI sessions.
%   res = buildFusiAnatomyNaive(opts, sessionvols, voxelsize_mm) builds a
%   within-mouse anatomy scan from many freehand-repositioned fUSI sessions the
%   simple, robust way: instead of the non-uniform super-resolution splat of
%   buildFusiAnatomy, EVERY session is upsampled once to the atlas resolution,
%   the rigid session-to-template transforms are estimated on that common grid,
%   and the aligned sessions are combined (median by default). It is the
%   first-class version of the "median(anatomy.sessionvolumes)" step that the
%   demo already relies on, hardened so a single bad session cannot abort or
%   poison the build:
%     1) UPSAMPLE: each session is resampled to isotropic atlas resolution
%        exactly once (no per-pass re-resampling) and equalized to a common
%        median so no session dominates the average.
%     2) SEED: the user picks a seed session (selectSeedSession, cached in
%        <opts.savepath>/seed_session_for_anatomy.txt, shared with
%        buildFusiAnatomy). The seed, on the atlas grid, anchors "seed space".
%     3) GROUPWISE REGISTRATION (Npass passes): every other session is rigidly
%        registered to the current template (the seed on pass 1) with elastix,
%        WARM-STARTED from the previous pass's transform on later passes. Each
%        fit is wrapped so it cannot crash the run, and only fits that pass an
%        alignment check (enough overlap with, and positive correlation to, the
%        template) are accepted; the rest fall back to their best-known
%        transform and are left out of the template. The template is rebuilt
%        each pass as the median over the sessions that passed, so one bad
%        session neither aborts the build nor corrupts the target the others
%        register to.
%     4) COMBINE: the final anatomy is the median (or mean) over the sessions
%        that passed, with the few uncovered voxels nearest-filled, written to
%        <opts.savepath>/<opts.mousename>_anatomy_naive.mat.
%
%   Like buildFusiAnatomy the template stays in the NATIVE session orientation
%   (only the sampling changes), so res.volume/res.voxelsize_mm feed straight
%   into setupFusiOptions.
%
%   Inputs:
%     opts         - struct with fields savepath, mousename and atlas_res.
%     sessionvols  - [d1 d2 d3 Nsessions] per-session volumes (median over time
%                    of each recording). Coronal slices run along dim 3.
%     voxelsize_mm - [d1 d2 d3] voxel size of sessionvols, in mm.
%     sessionnames - (optional) cellstr of Nsessions names.
%     savesessions - (optional, default true) also store every session warped
%                    into seed space (at atlas resolution) for QC.
%     Npass        - (optional, default 3) number of groupwise passes. Npass=1
%                    registers to the seed once with no refinement.
%     combinemode  - (optional, default 'median') 'median' or 'mean' central
%                    tendency used to build the template and final anatomy.
%                    'median' is robust to a residual mis-aligned session.
%     transformtype- (optional, default 'rigid') 'rigid' or 'affine' elastix
%                    model used to align each session. 'affine' also absorbs a
%                    global scale/shear between freehand sessions (e.g. a
%                    slightly wrong or session-varying voxelsize_mm) that a rigid
%                    fit cannot. A diverged affine fit that grossly rescales the
%                    brain is caught by the alignment check and dropped.
%
%   Output (also the contents of <mousename>_anatomy_naive.mat):
%     res.volume               - [D1 D2 D3] single anatomy at atlas resolution,
%                                NaN-free (holes filled), for setupFusiOptions.
%     res.voxelsize_mm         - resolution of res.volume (== opts.atlas_res)
%     res.seed_session         - index of the seed session
%     res.session_names        - names of all sessions, in order
%     res.session_voxelsize_mm - resolution of the input volumes
%     res.npass                - groupwise passes run
%     res.combine              - 'median' or 'mean' used
%     res.transformtype        - 'rigid' or 'affine' model used
%     res.tforms               - Nsessions cell of session->seed transforms
%     res.valid                - Nsessions logical: session met the alignment
%                                check and contributed to the anatomy
%     res.alignment_corr       - Nsessions final correlation with the template
%     res.session_coverage     - Nsessions final fraction of the template brain
%                                the warped session covered
%     res.sessionvolumes       - [D1 D2 D3 Nsessions] sessions warped into seed
%                                space at atlas resolution (QC only, optional)
%
%   See also BUILDFUSIANATOMY, SELECTSEEDSESSION, SETUPFUSIOPTIONS.

Nsessions = size(sessionvols, 4);

if nargin < 4 || isempty(sessionnames)
    sessionnames = arrayfun(@(x) sprintf('session%02d', x), 1:Nsessions, 'UniformOutput', false);
end
if nargin < 5 || isempty(savesessions), savesessions = true;     end
if nargin < 6 || isempty(Npass),        Npass        = 3;        end
if nargin < 7 || isempty(combinemode),  combinemode  = 'median'; end
if nargin < 8 || isempty(transformtype), transformtype = 'rigid'; end
combinemode   = validatestring(combinemode,   {'median', 'mean'});
transformtype = validatestring(transformtype, {'rigid', 'affine'});
voxelsize_mm = double(voxelsize_mm(:)).';
atlasres     = double(opts.atlas_res(:)).';
if isscalar(atlasres), atlasres = atlasres([1 1 1]); end

% Alignment acceptance thresholds. Data-driven, with no rotation/translation
% priors: a fitted session is kept only if, warped into seed space, it still
% covers most of the template brain and is positively correlated with it. A
% diverged or wrong-optimum elastix fit fails these and is dropped from the
% average instead of poisoning it (buildFusiAnatomy rebuilds its template from
% every session, so one bad fit propagates to all the others on the next pass).
mincorr = 0.30;   % Pearson corr with the template over their joint support
mincov  = 0.50;   % >= 50% of the template brain voxels covered by the session

makeNewDir(opts.savepath);
%==========================================================================
% (1) UPSAMPLE every session to isotropic atlas resolution ONCE. This is the
% "naive" step: no super-resolution grid, everything lives on the atlas grid,
% and the resampling is not repeated every pass (buildFusiAnatomy re-runs
% imresize3 on each session in every pass and again for QC).
facvol    = voxelsize_mm./atlasres;
sessatlas = cell(Nsessions, 1);
for isess = 1:Nsessions
    v = sanitizeVol(resampleSession(sessionvols(:, :, :, isess), facvol));
    m = median(v(v > 0), 'omitnan');          % isempty guards an all-zero session
    if isempty(m) || ~isfinite(m) || m <= 0, m = 1; end
    sessatlas{isess} = v./m;                  % equalized so none dominates
end
%==========================================================================
% (2) SEED (shares the cache file with buildFusiAnatomy)
iseed   = selectSeedSession(sessionvols, opts.savepath, sessionnames, voxelsize_mm);
seedvol = sessatlas{iseed};
Rseed   = imref3d(size(seedvol));
sz      = size(seedvol);
fprintf(['Building naive (atlas-resolution) anatomy for %s: %d sessions, ', ...
    'seed %d (%s), %d pass(es), %s combine, %s transform.\n'], ...
    opts.mousename, Nsessions, iseed, sessionnames{iseed}, Npass, combinemode, transformtype);
%==========================================================================
regpath = fullfile(opts.savepath, 'anatomy_registration_naive');
makeNewDir(regpath);

% warm-started passes are already roughly aligned, so elastix fits at a single
% full scale (no coarse pyramid that could jump away from the warm start)
residualparams = struct('NumberOfResolutions', 1, 'ImagePyramidSchedule', [1 1 1], ...
    'MaximumNumberOfIterations', 1000);

tforms    = repmat({makeTform(eye(4), transformtype)}, Nsessions, 1);
aligned   = nan([sz Nsessions], 'single');    % every session warped into seed space
validmask = false(Nsessions, 1);
corrs     = zeros(Nsessions, 1);
covs      = zeros(Nsessions, 1);
%==========================================================================
% (3) GROUPWISE PASSES: register to the evolving median template, keep only the
% fits that pass the alignment check, rebuild the template from those.
fixedfine = seedvol;                          % pass-1 registration/validation target
for pass = 1:Npass
    validmask(:) = false;                     % re-evaluated every pass
    for isess = 1:Nsessions
        %------------------------------------------------------------------
        if isess == iseed
            tforms{isess}         = makeTform(eye(4), transformtype);   % seed anchors the frame
            aligned(:, :, :, isess) = seedvol;
            validmask(isess)      = true;
            corrs(isess)          = 1;
            covs(isess)           = 1;
            continue;
        end
        %------------------------------------------------------------------
        outpath = fullfile(regpath, sessionnames{isess});
        makeNewDir(outpath);
        volatl  = sessatlas{isess};
        Ratl    = imref3d(size(volatl));
        %------------------------------------------------------------------
        % fit a candidate transform: cold to the seed on pass 1, warm-started
        % from the previous pass afterwards. A cold fit of a freehand session
        % (far from the template) can diverge, so it is wrapped: on any elastix
        % error we keep the previous transform and carry on instead of aborting
        % the whole multi-session build.
        Tnew      = [];
        cachefile = '';
        try
            if pass == 1
                [Tnew, cachefile] = registerSessionToTarget(fixedfine, volatl, atlasres(1), outpath, pass, transformtype);
            else
                Tprev             = tforms{isess};
                sesspre           = imwarp(volatl, Ratl, Tprev, 'OutputView', Rseed, 'FillValues', 0);
                [Tcorr, cachefile] = registerSessionToTarget(fixedfine, sesspre, atlasres(1), outpath, pass, transformtype, residualparams);
                Tnew              = makeTform(Tcorr.A * Tprev.A, transformtype);
            end
        catch ME
            warning('buildFusiAnatomyNaive:regFailed', ...
                'Session %d (%s) pass %d: elastix failed (%s). Keeping previous transform.', ...
                isess, sessionnames{isess}, pass, ME.message);
        end
        %------------------------------------------------------------------
        % never regress: score the incumbent transform and the new fit against
        % the current template and keep whichever is a better VALID alignment.
        [warpbest, okbest, ccbest, covbest] = scoreTform(volatl, Ratl, tforms{isess}, Rseed, fixedfine, mincorr, mincov);
        if ~isempty(Tnew)
            [warpnew, oknew, ccnew, covnew] = scoreTform(volatl, Ratl, Tnew, Rseed, fixedfine, mincorr, mincov);
            takenew = (oknew && (~okbest || ccnew > ccbest)) || (~oknew && ~okbest && ccnew > ccbest);
            if takenew
                tforms{isess} = Tnew;
                warpbest = warpnew; okbest = oknew; ccbest = ccnew; covbest = covnew;
            elseif ~oknew && ~isempty(cachefile) && exist(cachefile, 'file')
                % this pass's fit is not usable; drop its cache so a rerun refits
                % (elastix sampling is stochastic) rather than reloading the bad
                % fit. A valid-but-not-better fit is kept cached for reproducibility.
                delete(cachefile);
            end
        end
        %------------------------------------------------------------------
        aligned(:, :, :, isess) = warpbest;
        validmask(isess)        = okbest;
        corrs(isess)            = ccbest;
        covs(isess)             = covbest;
        if ~okbest
            warning('buildFusiAnatomyNaive:poorFit', ...
                'Session %d (%s) pass %d left out (corr %.2f, coverage %.0f%%).', ...
                isess, sessionnames{isess}, pass, ccbest, 100*covbest);
        end
        %------------------------------------------------------------------
    end
    template  = combineAligned(aligned, validmask, combinemode);
    fixedfine = sanitizeVol(fillMissingVolume(template));   % clean target for next pass
    fprintf('  Pass %d/%d: %d/%d sessions aligned (median corr %.2f).\n', ...
        pass, Npass, nnz(validmask), Nsessions, median(corrs(validmask)));
end
%==========================================================================
% (4) COMBINE + close the leftover holes
templatefull = fillMissingVolume(combineAligned(aligned, validmask, combinemode));
nrej         = Nsessions - nnz(validmask);
if nrej > 0
    warning('buildFusiAnatomyNaive:rejected', ...
        '%d of %d sessions did not meet the alignment threshold and are left out of the anatomy: %s.', ...
        nrej, Nsessions, strjoin(sessionnames(~validmask), ', '));
end

res = struct();
res.volume               = templatefull;
res.voxelsize_mm         = atlasres;
res.seed_session         = iseed;
res.session_names        = sessionnames(:);
res.session_voxelsize_mm = voxelsize_mm;
res.npass                = Npass;
res.combine              = combinemode;
res.transformtype        = transformtype;
res.tforms               = tforms;
res.valid                = validmask;
res.alignment_corr       = corrs;
res.session_coverage     = covs;
if savesessions
    res.sessionvolumes = single(aligned);
end
%==========================================================================
% QC: overlay each session, warped into seed space, on the final template
for isess = 1:Nsessions
    if isess == iseed, continue; end
    cf = visualizeTransformMatch(forDisplay(sqrt(templatefull)), forDisplay(sqrt(aligned(:, :, :, isess))), 3);
    print(cf, fullfile(regpath, sessionnames{isess}, 'rigid_registration_to_template'), '-dpng')
    close(cf);
end
%==========================================================================
savename = fullfile(opts.savepath, sprintf('%s_anatomy_naive.mat', opts.mousename));
save(savename, '-struct', 'res', '-v7.3');
fprintf('Naive anatomy volume [%s] at %.0f um saved to %s\n', ...
    num2str(size(res.volume)), atlasres(1)*1e3, savename);
%==========================================================================
% one overview figure of the anatomy against the seed alone
cf = visualizeTransformMatch(forDisplay(sqrt(seedvol)), forDisplay(sqrt(templatefull)), 3);
print(cf, fullfile(opts.savepath, sprintf('%s_anatomy_naive', opts.mousename)), '-dpng')
close(cf);
%==========================================================================
end

function [tform, targetfile] = registerSessionToTarget(targetvol, sessvol, volscale, outpath, pass, transformtype, paramoverrides)
%REGISTERSESSIONTOTARGET Elastix rigid/affine fit of a session to the template.
%   Registers the target (seed on pass 1, the evolving template afterwards)
%   onto the session and inverts, so the returned transform maps the session
%   into seed space. transformtype is 'rigid' (rigidtform3d) or 'affine'
%   (affinetform3d); the two write to separate caches so they never collide.
%   Cached per session and per pass so reruns are cheap; the cache path is
%   returned so the caller can invalidate a fit that fails its alignment check.
%   paramoverrides (optional) forwards elastix overrides, e.g. a single
%   full-scale resolution for the warm-started passes.
if nargin < 7, paramoverrides = struct(); end
targetfile = fullfile(outpath, sprintf('%s_to_template_pass%d.txt', transformtype, pass));
if exist(targetfile, 'file')
    fprintf('    Found %s transform (pass %d) in %s, loading...\n', transformtype, pass, outpath)
else
    fprintf('    Fitting %s transform (pass %d)...\n', transformtype, pass)
    if strcmp(transformtype, 'affine')
        [~, ~, tformpath, ~] = performElastixAffineRegistration(targetvol, sessvol, volscale, outpath, paramoverrides);
    else
        [~, ~, tformpath, ~] = performElastixRigidRegistration(targetvol, sessvol, volscale, outpath, paramoverrides);
    end
    movefile(tformpath, targetfile);
end
tform = parse_elastix_tform_all(targetfile);
tform = makeTform(tform.A, transformtype);
tform = tform.invert();
end

function tf = makeTform(A, transformtype)
%MAKETFORM Wrap a 4x4 matrix as the chosen geometric transform object.
%   'affine' keeps the full 12-DOF matrix; 'rigid' forces it back onto the
%   rigid (rotation+translation) manifold, dropping numeric non-orthogonality.
if strcmp(transformtype, 'affine')
    tf = affinetform3d(A);
else
    tf = rigidtform3d(A);
end
end

function [warpvol, ok, cc, covfrac] = scoreTform(vol, Rin, tform, Rseed, target, mincorr, mincov)
%SCORETFORM Warp a session into seed space and score it against the template.
%   Returns the warped volume plus whether it passes the alignment check, its
%   correlation with the template and the fraction of the template brain it
%   covers. Used to accept/reject each fit and to pick the better of two.
warpvol            = imwarp(vol, Rin, tform, 'OutputView', Rseed, 'FillValues', nan);
warpvol(warpvol < 0) = 0;
[ok, cc, covfrac]  = checkAlignment(warpvol, target, mincorr, mincov);
% also reject a transform that grossly rescales the volume: a diverged affine
% fit collapses or explodes the brain (det far from 1). Rigid fits have det ~ 1
% and always pass this, so the guard only bites in 'affine' mode.
dscale = abs(det(tform.A(1:3, 1:3)));
ok     = ok && dscale > 1/3 && dscale < 3;
end

function [ok, cc, covfrac] = checkAlignment(warpvol, target, mincorr, mincov)
%CHECKALIGNMENT Overlap fraction and Pearson correlation over the joint support.
%   A well-aligned session covers most of the template brain (target > 0) with
%   real signal and is positively correlated with it there. A diverged fit maps
%   the brain out of view (low coverage) or onto unrelated structure (low/neg
%   correlation), so it fails one of the two thresholds.
brain   = nnz(target > 0);
support = isfinite(warpvol) & isfinite(target) & target > 0;
covfrac = nnz(support & warpvol > 0)/max(brain, 1);   % scalar
if nnz(support) < 100
    cc = 0; ok = false; return;
end
a  = double(warpvol(support)); a = a - mean(a);
b  = double(target(support));  b = b - mean(b);
den = norm(a)*norm(b);                                % scalars
if den <= 0
    cc = 0;
else
    cc = (a.'*b)/den;
end
ok = isfinite(cc) && cc >= mincorr && covfrac >= mincov;
end

function vol = combineAligned(aligned, mask, mode)
%COMBINEALIGNED Median/mean over the sessions that passed, NaN-omitting.
idx = find(mask);
if isempty(idx)
    error('buildFusiAnatomyNaive:noValidSessions', ...
        'No sessions passed the alignment check; cannot build a template.');
end
sub = aligned(:, :, :, idx);
if strcmp(mode, 'median')
    vol = median(sub, 4, 'omitnan');
else
    vol = mean(sub, 4, 'omitnan');
end
end

function vol = sanitizeVol(vol)
%SANITIZEVOL Force a volume to finite, non-negative single precision.
%   Real I./Inorm data can hold Inf/NaN (where Inorm == 0) and small negatives;
%   either makes elastix diverge or crash, so clean them before registration.
vol = single(vol);
vol(~isfinite(vol)) = 0;
vol(vol < 0)        = 0;
end

function vol = resampleSession(vol, facvol)
%RESAMPLESESSION Bring one session to atlas sampling, as single.
%   Intensity is not normalized here (elastix uses scale-invariant mutual
%   information); the per-session median equalization is applied by the caller.
vol = imresize3(single(vol), 'Scale', facvol);
end

function vol = fillMissingVolume(vol)
%FILLMISSINGVOLUME 3-D nearest-neighbour fill of NaN voxels.
%   The N-D analogue of fillmissing(...,'nearest'): every NaN voxel takes the
%   value of the closest voxel that had data. Guarded so an all-NaN volume (no
%   coverage at all) is returned unchanged instead of raising a bwdist index-0
%   error, as buildFusiAnatomy's version would.
nanmask = isnan(vol);
if any(nanmask, 'all') && ~all(nanmask, 'all')
    [~, idx]     = bwdist(~nanmask);
    vol(nanmask) = vol(idx(nanmask));
end
end

function im = forDisplay(vol)
%FORDISPLAY uint8 copy scaled to its 1-99% range, for the imfuse overlays in
%   visualizeTransformMatch. quantile ignores the NaN FOV border, and NaNs map
%   to 0 (black) in the uint8 cast, so no extra masking is needed.
vol    = single(vol);
imlims = quantile(vol, [0.01 0.99], 'all');
im     = uint8(255*(vol - imlims(1))/range(imlims));
end

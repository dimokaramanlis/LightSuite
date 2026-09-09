% LS_ANALYZE_FUSI  End-to-end functional ultrasound (fUS) workflow.
%
% This demo takes a set of repeated, freehand-repositioned fUS scans of one
% mouse and ends with functional maps in Allen CCF space (and, optionally, on
% the Allen cortical flatmap). It is written to be adapted: every path is a
% placeholder, and the only assumption about your data is that you can produce
%   * a set of per-session anatomical volumes (one 3-D volume per session), and
%   * the functional recordings themselves (a 4-D volume, time last).
%
% The pipeline has five stages:
%   1. ACROSS-SESSION RIGID ALIGNMENT. Every session is upsampled to the atlas
%      resolution and rigidly aligned to a user-picked seed session. The mean
%      of the aligned sessions is the mouse's "anatomy scan" and defines the
%      common seed space.                            -> buildFusiAnatomy
%   2. ORIENTATION. You tell LightSuite how the sample axes map onto the atlas
%      axes (once per mouse, cached).                -> setupFusiOptions
%   3. LANDMARKS. You place matched control points on the anatomy and on the
%      atlas in a side-by-side GUI.                  -> matchControlPoints_minimal
%   4. THE FIT. An affine + B-spline registration is optimized against BOTH the
%      image similarity and your landmarks.          -> multiobjRegistrationFusi
%   5. APPLYING THE TRANSFORM. A functional recording is rigidly aligned to the
%      anatomy, and that rigid step chained with the anatomy->atlas warp brings
%      any map or timeseries into atlas space.
%                                 -> fusiRecordingToAnatomy + applyFusiTransforms
%
% Stages 1-4 are run ONCE per mouse; their results are cached on disk, so
% re-running the script skips straight past the GUIs. Stage 5 is run once per
% recording.
%
% Prerequisites: elastix on the system PATH, the Allen CCF on the MATLAB path,
% and the fUS vascular atlas (see src/fusi/prepare_fusi_atlas.m). See the
% "Functional ultrasound" page of the documentation for the full description.

%==========================================================================
%% (setup) atlas, paths and options
%==========================================================================
% The fUS workflow registers against a VASCULAR atlas: the Allen CCF geometry
% carrying a vascular-contrast template, so that image similarity is computed
% between two images of the same modality (vessels against vessels).
[tvvessel, avvessel, parcelinfo, atlas_res] = loadAtlasInfo('allen2020fusi_50um');

opts            = struct();
opts.atlas      = tvvessel;      % template volume (vascular contrast)
opts.annotation = avvessel;      % annotation / label volume
opts.atlas_res  = atlas_res;     % [1 1 1]*0.05 mm
opts.parcelinfo = parcelinfo;    % Allen parcellation table

pathrawdata     = 'C:\data\fusi';          % EDIT ME: root of your raw data
opts.mousename  = 'mouse01';               % EDIT ME
opts.savepath   = fullfile(pathrawdata, opts.mousename, 'lightsuite');
makeNewDir(opts.savepath);

% Voxel size of the RAW sessions, in mm, along the raw session dimensions.
% Coronal slices are expected to run along dimension 3.
pxsizesession   = [0.1971, 0.15, 0.15];    % EDIT ME
finshape        = [36, 64, 54];            % EDIT ME: [d1 d2 d3] of one session

%==========================================================================
%% (1a) load one anatomical volume per session
%==========================================================================
% Each session contributes a single 3-D volume: a robust average over time of
% that recording (a "session template"). How you get it depends on your file
% format, so replace the loop below with your own loader.
%
% For recordings stored as a *_FUS.mat file holding I = [nvox x nframes],
% loadFusiSessionTemplate does this for you (and caches the result next to the
% raw file, so the second run is instant).
allpaths     = dir(fullfile(pathrawdata, opts.mousename, '**', '*FUS.mat'));
allsessions  = fullfile({allpaths(:).folder}', {allpaths(:).name}');
sessionnames = erase({allpaths(:).name}', '_FUS.mat');
Nsessions    = numel(allsessions);

sessionvols  = nan([finshape, Nsessions], 'single');
tic; msg = [];
for isess = 1:Nsessions
    sessionvols(:, :, :, isess) = loadFusiSessionTemplate(allsessions{isess}, finshape);
    fprintf(repmat('\b', 1, numel(msg)));
    msg = sprintf('Session %d/%d done. Took %2.2f s\n', isess, Nsessions, toc);
    fprintf(msg);
end

%==========================================================================
%% (1b) (manual on first run) rigid alignment across sessions -> anatomy scan
%==========================================================================
% buildFusiAnatomy upsamples every session to the atlas resolution, asks you to
% pick a SEED session (a GUI showing all sessions side by side; the choice is
% cached in seed_session_for_anatomy.txt), rigidly registers every other
% session to that seed with elastix, and averages the aligned volumes.
%
% Pick a seed that is well centred, artefact-free and covers as much of the
% brain as possible: it defines the common "seed space" that everything else is
% expressed in, so a bad seed costs you on every later step.
%
% Per-session transforms are cached under <savepath>/anatomy_registration/, so
% the elastix fits only run the first time.
anatomy = buildFusiAnatomy(opts, sessionvols, pxsizesession, sessionnames);

% anatomy.volume         - the averaged anatomy, at atlas resolution
% anatomy.tforms{i}      - session i -> seed rigid transform (reused in stage 5)
% anatomy.seed_session   - index of the seed session
% anatomy.sessionvolumes - every session warped into seed space (QC)
%
% Averaging is the robust default, but it blurs vessels when the alignment is
% imperfect. The seed session on its own is sharper and is a good alternative:
%   volanatomy = anatomy.sessionvolumes(:, :, :, anatomy.seed_session);
volanatomy = anatomy.volume;

% QC: scroll through the aligned sessions and check that they superimpose.
% for ii = 1:size(anatomy.sessionvolumes, 3)
%     imagesc(sqrt(squeeze(anatomy.sessionvolumes(:, :, ii, :)))); pause;
% end

%==========================================================================
%% (2) (manual on first run) brain orientation
%==========================================================================
% setupFusiOptions asks (via a GUI) how the sample axes map onto the atlas
% axes, stores the permutation in brain_orientation.txt, attaches the anatomy
% to opts as the registration sample, and writes <savepath>/regopts.mat.
% On later runs the stored orientation is reused and no GUI appears.
%
% The anatomy from stage 1 is already AT the atlas resolution, so the voxel
% size to pass here is opts.atlas_res -- not pxsizesession.
opts = setupFusiOptions(volanatomy, opts.atlas_res, opts);

%==========================================================================
%% (3) (manual) place matched control points
%==========================================================================
% Side-by-side GUI: your anatomy on the left, the atlas on the right, both
% sliced along dimension 1. Click a landmark on one side, then its counterpart
% on the other. Vessel bifurcations, the midline, the brain outline and the
% ventricles all make good landmarks. Points are saved to
% <savepath>/control_points_minimal.mat and reloaded when you reopen the GUI.
%
% At least five pairs spanning all three dimensions are required. In the paper a
% median of 219 landmarks per mouse was used; spreading them over the full
% anteroposterior extent matters more than the raw count.
matchControlPoints_minimal(opts);

%==========================================================================
%% (4) (auto) optimize the anatomy -> atlas transform
%==========================================================================
% multiobjRegistrationFusi fits an affine and then a B-spline transform,
% optimizing a weighted sum of image similarity (mutual information between the
% vascular atlas and your anatomy) and the distance between your matched
% control points. It writes <savepath>/transform_params.mat, plus per-dimension
% overlay PNGs to check the result.
wtpoints                   = 0.1;  % weight of the landmark term vs image similarity
opts.bspline_spatial_scale = 1.6;  % mm; smaller = more local deformation
opts.n_histogram_bins      = 48;   % bins for the mutual-information estimate
multiobjRegistrationFusi(opts, wtpoints, false);

% If the fit drifts away from your landmarks, raise wtpoints (0.1 -> 1 -> 2).
% If it looks over-warped, raise bspline_spatial_scale. Always inspect the
% <mouse>_dim*_affine_registration.png and *_bspline_registration.png overlays
% before trusting anything downstream.

%==========================================================================
%% (5a) load a functional recording
%==========================================================================
% Any 4-D array [d1 d2 d3 nframes] in the raw session geometry works. Here we
% take one recording of this mouse; adapt the loading to your own format.
irec     = 1;
recname  = sessionnames{irec};
recdata  = load(allsessions{irec});                      % EDIT ME
scandata = reshape(single(recdata.I), [finshape, []]);   % [d1 d2 d3 nframes]

%==========================================================================
%% (5b) (auto) rigidly bring the recording into anatomy (seed) space
%==========================================================================
% A functional recording sits wherever the probe was on that day, so it first
% has to be rigidly aligned to the anatomy. fusiRecordingToAnatomy fits that
% rigid transform against the recording's own time-median, caches it under
% <savepath>/<recname>/, and prints a QC overlay.
tformfuntoanatomy = fusiRecordingToAnatomy(opts, scandata, pxsizesession, recname);

% If the recording is one of the sessions that built the anatomy, its rigid
% transform is ALREADY known and does not need re-fitting -- reuse it:
%   tformfuntoanatomy = struct('tformrigid', anatomy.tforms{irec}, ...
%                              'Ranatomy',   imref3d(size(volanatomy)));

%==========================================================================
%% (5c) (auto) carry a functional MAP into atlas space
%==========================================================================
% Compute whatever map you like in native session space -- a stimulus
% correlation map, a GLM beta or T-score, a difference between conditions.
% Here: correlation with a boxcar stimulus between frames t0 and t1, convolved
% with a hemodynamic response function.
t0 = 30; t1 = 40;
scanfus = struct('Data', scandata, 'VoxelSize', pxsizesession*1e3);
outmap  = mapCorrelation(scanfus, t0, t1);

% applyFusiTransforms chains the rigid step (recording -> anatomy) with the
% B-spline + affine warp (anatomy -> atlas). savefullvols = true also returns
% the warped volumes; set it to false when the per-area averages are enough
% (much less memory).
rescorr = applyFusiTransforms(opts, outmap.Data, pxsizesession, true, tformfuntoanatomy);

% rescorr.vreg        - the map in atlas space, full non-rigid warp
% rescorr.vregaff     - the same map, affine only (a useful sanity check)
% rescorr.areasignals - [nArea x 1] area averages of the warped map
% rescorr.groupidx    - Allen parcellation ids labelling those rows
% rescorr.areavols    - volume of each area, in mm^3
save(fullfile(opts.savepath, recname, 'corr_registered.mat'), '-struct', 'rescorr', '-v7.3');

% quick look: coronal montage of the atlas-space map over the atlas template
cf = plotFusiActivationMontage(rescorr.vreg, struct( ...
    'underlay',  single(opts.atlas), ...
    'thresh',    0.1, ...
    'sliceaxis', 1, ...
    'sliceidx',  round(linspace(1, size(rescorr.vreg, 1), 24)), ...
    'title',     sprintf('%s - %s stimulus correlation', opts.mousename, recname)));

%==========================================================================
%% (5d) (auto) carry a whole TIMESERIES into atlas space
%==========================================================================
% The same call works on a 4-D array: every frame is warped with the same
% (precomputed) deformation field. Full atlas-space movies are large, so keep
% savefullvols = false unless you really need the voxelwise timeseries -- the
% per-area timecourses in .areasignals are usually what you are after.
restime = applyFusiTransforms(opts, scandata, pxsizesession, false, tformfuntoanatomy);

% restime.areasignals - [nArea x nframes] atlas-area timecourses
save(fullfile(opts.savepath, recname, 'time_registered.mat'), '-struct', 'restime', '-v7.3');

% Example: the timecourse of one named structure. The rows of .areasignals are
% parcellation LEAVES, labelled by .groupidx; a structure such as VISp covers
% several of them (its layers), so pool them weighted by their volume.
%
% parcelinfo lists every leaf once per term set ('substructure', 'structure',
% 'division', ...), so pick the level you want before matching an acronym.
isstruct = strcmp(opts.parcelinfo.parcellation_term_set_name, 'structure');
leafidx  = opts.parcelinfo.parcellation_index( ...
    isstruct & strcmpi(opts.parcelinfo.parcellation_term_acronym, 'VISp'));
irows    = ismember(restime.groupidx, leafidx);
if any(irows)
    w  = restime.areavols(irows);
    tc = sum(w .* restime.areasignals(irows, :), 1, 'omitnan') / sum(w);
    figure; plot(tc); xlabel('frame'); ylabel('signal'); title('VISp');
end

%==========================================================================
%% (6) (optional) project an atlas-space volume onto the cortical flatmap
%==========================================================================
% Any volume that is already in CCF space (AP x DV x ML, isotropic) can be
% flattened onto the Allen cortical surface: for every point of the flatmap the
% volume is averaged along the corresponding cortical streamline.
%
% The projection itself runs in PYTHON (allensdk + ccf_streamlines) via
% src/fusi/py/fusi_flatmap_project.py; MATLAB only hands over the volume and
% reads the 2-D result back, so the interpreter can live in any conda/venv
% environment. You also need the flatmap resources (flatmap_butterfly.h5/.nrrd,
% surface_paths_10_v3.h5, labelDescription_ITKSNAPColor.txt, manifest.json).
fmopts = struct();
fmopts.pythonexe   = 'C:\miniconda3\envs\allensdk\python.exe';  % EDIT ME
fmopts.resourcedir = 'C:\AllenAtlas\flatmaps';                  % EDIT ME
fmopts.inputres_um = 50;        % our atlas-space volumes are 50 um isotropic
fmopts.hemisphere  = 'left';
fmopts.kind        = 'mean';    % streamline average ('max' is also available)

vol = rescorr.vreg;

% Optional: symmetrize across the midline before projecting. fUS coverage is
% rarely identical on the two sides, and mirroring fills one hemisphere with
% whatever the other one saw. Skip this if left/right differences matter.
Nmid = size(vol, 3)/2;
vol  = (vol(:, :, 1:Nmid) + flip(vol(:, :, Nmid+1:end), 3))/2;
vol  = cat(3, vol, flip(vol, 3));

fm = fusiVolumeToFlatmap(vol, fmopts);

% fm.flatmap          - 2-D projection
% fm.boundaries       - area-boundary overlay image (same size), or []
% fm.regionBoundaries - per-area outlines (.names, .coords), to draw a subset
figure('Color', 'w');
ax   = axes();
cmax = max(abs(fm.flatmap), [], 'all', 'omitnan')*0.8;
imagesc(ax, fm.flatmap.', [-1 1]*cmax);
axis(ax, 'image', 'off'); colormap(ax, turbo); colorbar(ax);
hold(ax, 'on');
if ~isempty(fm.boundaries)
    [yy, xx] = find(fm.boundaries.');
    plot(ax, xx, yy, '.k', 'MarkerSize', 1);
end
title(ax, sprintf('%s - flatmap projection', opts.mousename));

function atl = warpFusiSessionToAtlas(optA, tform, voxelsize_mm, sr, opts)
%WARPFUSISESSIONTOATLAS Move a session's activation maps/timecourses to atlas.
%   atl = warpFusiSessionToAtlas(optA, tform, voxelsize_mm, sr) brings the
%   native-space results in sr (from fusiSessionActivationMaps) into Allen
%   atlas space, reusing the registration already computed by the pipeline:
%     1) the session->seed rigid transform is taken from the mouse anatomy
%        (<mouse>_anatomy.mat .tforms{i}, built by buildFusiAnatomy) - it is
%        NOT re-fitted. The maps are upsampled to atlas resolution and warped
%        into seed / anatomy space with it (same step as buildFusiAnatomy).
%     2) applyFusiTransforms then performs only the seed->atlas B-spline +
%        affine warp (data is already in anatomy space, so no rigid step), and
%        returns BOTH the full non-rigid and the affine-only result.
%   Maps are passed as the frames of a 4-D volume, so one warm deformation
%   field serves them all.
%
%   Input:
%     optA          atlas opts from regopts.mat (atlas, annotation, atlas_res,
%                   parcelinfo, sample [seed-space anatomy], sample_res,
%                   permute_sample_to_atlas, savepath+transform_params.mat).
%     tform         session->seed rigid transform (rigidtform3d), i.e.
%                   <mouse>_anatomy.mat .tforms{i} for this session.
%     voxelsize_mm  session voxel size [d1 d2 d3] in mm.
%     sr            session result from fusiSessionActivationMaps.
%     opts          (optional): .finshape (default [36 64 54]),
%                   .areatimecourses (default true), .verbose (default true).
%
%   Output struct atl (each non-rigid field has an affine-only twin, so the
%   two registrations can be compared):
%     .maps/.mapsAff        struct of atlas-space volumes [atlassize]:
%                  combined_corr, object_corr, scrambled_corr, combined_tscore.
%                  .maps = full non-rigid, .mapsAff = affine only.
%     .mapAreaVals/.mapAreaValsAff [nArea x 4] area-mean of the four maps.
%     .groupidx    Allen area ids matching the rows of mapAreaVals / areaTC
%     .areaTC/.areaTCaff [nArea x nlags x 3] per-area peristimulus timecourse
%                  for combined/object/scrambled, or [] if not requested.
%     .lags        [nlags x 1] (from sr.timecourses), or []
%     .mapnames    {'combined_corr','object_corr','scrambled_corr','combined_tscore'}
%
%   See also BUILDFUSIANATOMY, APPLYFUSITRANSFORMS, FUSISESSIONACTIVATIONMAPS.

if nargin < 5, opts = struct(); end
if ~isfield(opts, 'finshape') || isempty(opts.finshape), opts.finshape = [36 64 54]; end
if ~isfield(opts, 'areatimecourses') || isempty(opts.areatimecourses), opts.areatimecourses = true; end
if ~isfield(opts, 'verbose'), opts.verbose = true; end

Rseed  = imref3d(size(optA.sample));          % seed / anatomy grid (atlas res)
facvol = voxelsize_mm ./ optA.atlas_res;      % native -> atlas resolution

% (1) maps -> seed/anatomy space with the anatomy rigid transform
mapnames = {'combined_corr','object_corr','scrambled_corr','combined_tscore'};
mapstack = zeroNaN(cat(4, sr.maps.combined.corr, sr.maps.object.corr, ...
                          sr.maps.scrambled.corr, sr.maps.combined.tscore));
mapseed  = stackToSeed(mapstack, facvol, tform, Rseed);

% (2) seed -> atlas (B-spline + affine, both non-rigid and affine-only)
resM = applyNoRigid(optA, mapseed, true);

atl = struct();
atl.mapnames = mapnames;
atl.maps = struct(); atl.mapsAff = struct();
for k = 1:numel(mapnames)
    atl.maps.(mapnames{k})    = resM.vreg(:, :, :, k);
    atl.mapsAff.(mapnames{k}) = resM.vregaff(:, :, :, k);
end
atl.mapAreaVals    = resM.areasignals;      % [nArea x 4] non-rigid
atl.mapAreaValsAff = resM.areasignalsaff;   % [nArea x 4] affine only
atl.groupidx  = resM.groupidx;
atl.areaTC    = [];
atl.areaTCaff = [];
atl.lags      = [];

% (3) per-area peristimulus timecourses (optional)
if opts.areatimecourses && ~isempty(sr.timecourses)
    tc    = sr.timecourses;
    nlags = numel(tc.lags);
    vol   = @(m) reshape(m, [opts.finshape nlags]);   % [nvox x nlags] -> 4-D
    peristack = zeroNaN(cat(4, vol(tc.perVox.combined), ...
                               vol(tc.perVox.object), vol(tc.perVox.scrambled)));
    periseed  = stackToSeed(peristack, facvol, tform, Rseed);
    resT = applyNoRigid(optA, periseed, false);
    atl.areaTC    = reshape(resT.areasignals,    [], nlags, 3);   % non-rigid
    atl.areaTCaff = reshape(resT.areasignalsaff, [], nlags, 3);   % affine only
    atl.lags   = tc.lags;
end
if opts.verbose
    fprintf('  warped session to atlas (%d areas%s).\n', numel(atl.groupidx), ...
        ternary(isempty(atl.areaTC), '', ', + peristimulus timecourses'));
end
end

% -------------------------------------------------------------------------
function S = stackToSeed(V, facvol, tform, Rseed)
%STACKTOSEED Upsample each frame to atlas res and warp it into seed space.
%   Mirrors the session->seed step of buildFusiAnatomy, frame by frame.
K = size(V, 4);
S = nan([Rseed.ImageSize K], 'single');
for k = 1:K
    vk = imresize3(single(V(:, :, :, k)), 'Scale', facvol);
    S(:, :, :, k) = imwarp(vk, imref3d(size(vk)), tform, ...
        'OutputView', Rseed, 'FillValues', nan);
end
end

% -------------------------------------------------------------------------
function res = applyNoRigid(optA, seeddata, savevols)
%APPLYNORIGID seed->atlas warp; data is already in anatomy space (no rigid).
w = warning('off', 'applyFusiTransforms:noRigidTransform');
cleaner = onCleanup(@() warning(w));
res = applyFusiTransforms(optA, seeddata, optA.atlas_res, savevols);
end

% -------------------------------------------------------------------------
function v = zeroNaN(v)
v(~isfinite(v)) = 0;   % undefined (background/constant) voxels -> no signal
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end

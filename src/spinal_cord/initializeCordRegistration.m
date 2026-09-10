function opts = initializeCordRegistration(inputpath, varargin)
%INITIALIZECORDREGISTRATION Second half of the spinal cord initialisation.
%
%   OPTS = INITIALIZECORDREGISTRATION(INPUTPATH) straightens the cord using the
%   centre line agreed on in SPINAL_CORD_ALIGNER and brings the atlas onto it
%   with an affine transform, which is the state the control point GUI and the
%   B-spline step expect. INPUTPATH is the folder holding regopts.mat.
%
%   The counterpart for a brain is INITIALIZEREGISTRATION. It only has to find
%   the orientation of the sample and a similarity transform, because a brain
%   does not change shape. A cord does, so this runs in two stages:
%
%     1. Straightening. Every slice is translated and rotated so the section
%        centre sits at the middle of the atlas frame and the dorsoventral axis
%        points the same way throughout. What comes out is a cord with the same
%        cross-sectional frame as the atlas, still at its own length.
%     2. Atlas fitting. The atlas is stretched along the long axis to cover the
%        straightened sample, then refined with an image-based affine
%        registration (elastix).
%
%   Both stages are written to disk: the straightened volume as
%   'cord_straight_register_<res>um.tif' (OPTS.straightvolpath), the transforms
%   into regopts.mat, and two annotated PNGs showing the atlas on the sample
%   before and after the affine refinement.
%
%   Optional name-value arguments
%     'zcoverage'        - fraction of the straightened sample the atlas is
%                          stretched to cover initially (default 0.98).
%     'orientationdeg'   - orientation the dorsoventral axis is rotated to
%                          (default 90, anterior towards the top of the image).
%     'skipaffine'       - keep the initial stretch and skip the elastix affine
%                          refinement (default false).
%
%   See also PREPARECORDSAMPLEFORREGISTRATION, SPINAL_CORD_ALIGNER,
%   MATCHCONTROLPOINTS_UNIFIED, MULTIOBJCORDREGISTRATION.

%==========================================================================
p = inputParser;
addRequired(p,  'inputpath', @(x) ischar(x) || isstring(x) || isstruct(x));
addParameter(p, 'zcoverage',      0.98, @(x) isscalar(x) && x > 0 && x <= 1);
addParameter(p, 'orientationdeg',   90, @isscalar);
addParameter(p, 'skipaffine',    false, @(x) islogical(x) || isscalar(x));
parse(p, inputpath, varargin{:});
params = p.Results;
%==========================================================================
opts = loadRegOpts(inputpath);
assert(isfield(opts, 'cordvolpath'), 'initializeCordRegistration:notPrepared', ...
    'Run prepareCordSampleForRegistration before initializeCordRegistration.');
%==========================================================================
% the centre line the user signed off on
alignpath = fullfile(opts.savepath, 'spinal_alignment_opt.mat');
if ~exist(alignpath, 'file')
    error('initializeCordRegistration:noAlignment', ...
        ['No alignment found in %s. Run spinal_cord_aligner and press ''s'' to ' ...
         'save before initialising the registration.'], opts.savepath);
end
aligndata = load(alignpath);
align_out = aligndata.align_out;
%==========================================================================
fprintf('Loading the cord atlas at %d um... ', opts.registres); tic;
[tv, av, atlasinfo] = loadCordAtlasVolumes(opts);
fprintf('Done! Took %2.2f s.\n', toc);
opts.atlasres        = atlasinfo.atlasres;
opts.atlassize       = atlasinfo.regsize;
opts.atlassizenative = atlasinfo.nativesize;
%==========================================================================
% 1. straighten the sample
targetcent    = 0.5 * atlasinfo.regsize([2 1]);   % [x y] centre of the atlas frame
opts.straighten = cordStraightenParams(align_out, targetcent, params.orientationdeg);
tforms          = computeStraighteningTransforms(opts.straighten);
opts.slicetforms = tforms;

cordvol = readDownStack(opts.cordvolpath);
Nslices = size(cordvol, 3);
assert(Nslices == opts.straighten.Nslices, ...
    'initializeCordRegistration:sliceMismatch', ...
    ['The alignment covers %d slices but %s has %d. Re-run spinal_cord_aligner ' ...
     'on the current volume.'], opts.straighten.Nslices, opts.cordvolpath, Nslices);

fprintf('Applying straightening transforms... '); tic;
raout       = imref2d(atlasinfo.regsize([1 2]));
straightvol = transformCordImageSlices(cordvol, tforms, raout);
fprintf('Done! Took %2.2f s.\n', toc);
clear cordvol;
%==========================================================================
% 2. bring the atlas onto the straightened sample
% first the obvious part: the atlas is a whole cord, the sample is a piece of
% one, so scale and shift along the long axis to cover it
z_scale = Nslices * params.zcoverage / atlasinfo.regsize(3);
z_trans = Nslices/2 - z_scale * atlasinfo.regsize(3)/2;

T      = eye(4);
T(3,3) = z_scale;
T(3,4) = z_trans;
transinit = affinetform3d(T);

refsample = imref3d(size(straightvol));
refatlas  = imref3d(atlasinfo.regsize);

volmax  = single(quantile(straightvol, 0.999, 'all'));
volplot = uint8(255 * single(straightvol) / volmax);

avsim = imwarp(av, refatlas, transinit, 'nearest', 'OutputView', refsample);
cf    = plotCordAnnotation(volplot, avsim);
print(cf, fullfile(opts.savepath, 'registration_initial_similarity'), '-dpng');
close(cf);
%==========================================================================
% then refine it with the images themselves
transaff = transinit;
if ~params.skipaffine
    tvtemp   = medfilt3(tv);
    atlasuse = imwarp(tvtemp, refatlas, transinit, 'OutputView', refsample);

    [~, ~, tformpath, ~] = performElastixAffineRegistration(...
        atlasuse, straightvol, 1, opts.savepath);
    newtrans = affinetform3d(parse_elastix_tform(tformpath));
    transaff = affinetform3d(transinit.A * newtrans.A);

    avshow = imwarp(av, refatlas, transaff, 'nearest', 'OutputView', refsample);
    cf     = plotCordAnnotation(volplot, avshow);
    print(cf, fullfile(opts.savepath, 'registration_initial_affine'), '-dpng');
    close(cf);
end
%==========================================================================
% save the straightened volume next to the other volumes
opts.straightvolpath = fullfile(opts.savepath, ...
    sprintf('cord_straight_register_%dum.tif', opts.registres));
saveopts = struct('compress', 'lzw', 'message', false);
if exist(opts.straightvolpath, 'file')
    delete(opts.straightvolpath);
end
saveastiff(straightvol, opts.straightvolpath, saveopts);

opts.affine_atlas_to_samp = transaff;
opts.straightvolsize      = size(straightvol);
%==========================================================================
saveRegOpts(opts);
fprintf(['Registration initialised. Add control points with ' ...
    'matchControlPoints_unified, then run multiobjCordRegistration.\n']);
%==========================================================================
end

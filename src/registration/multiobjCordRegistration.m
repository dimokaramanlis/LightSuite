function params_pts_to_atlas = multiobjCordRegistration(inputpath, control_point_wt, varargin)
%MULTIOBJCORDREGISTRATION Non-linear registration of a spinal cord to the atlas.
%
%   TRSTRUCT = MULTIOBJCORDREGISTRATION(INPUTPATH, CONTROL_POINT_WT) runs the
%   B-spline step for a cord, the counterpart of MULTIOBJREGISTRATION for a
%   brain. INPUTPATH is the folder holding regopts.mat. CONTROL_POINT_WT is how
%   much the user's control points count against the image similarity term; 0
%   ignores them entirely.
%
%   It expects the sample to have been straightened and the atlas affinely
%   fitted to it by INITIALIZECORDREGISTRATION. With enough control points
%   (more than 16) the affine transform is refitted from them first, then
%   elastix optimises a B-spline against both the images and the points, and the
%   result is inverted so that sample coordinates can be pushed into the atlas.
%
%   TRSTRUCT = MULTIOBJCORDREGISTRATION(..., 'bspline_spatial_scale', S) sets
%   the B-spline grid spacing in mm (default 0.96): smaller values allow more
%   local bending but follow noise more readily.
%
%   The returned struct is also saved as 'transform_params.mat' in the sample
%   folder. It holds the whole chain from raw sample pixels to atlas voxels -
%   downsampling, permutation, cropping, straightening, B-spline, affine and the
%   rostrocaudal flip - which is what GENERATEREGISTEREDCORDVOLUME and
%   CORDPOINTSTOATLAS replay.
%
%   See also INITIALIZECORDREGISTRATION, GENERATEREGISTEREDCORDVOLUME,
%   TRANSFORMCORDPOINTSTOATLAS, MULTIOBJREGISTRATION.

%==========================================================================
p = inputParser;
addRequired(p,  'inputpath', @(x) ischar(x) || isstring(x) || isstruct(x));
addRequired(p,  'control_point_wt', @isscalar);
addParameter(p, 'bspline_spatial_scale', 0.96, @(x) isscalar(x) && x > 0);
parse(p, inputpath, control_point_wt, varargin{:});
params = p.Results;
%==========================================================================
fprintf('Loading data volume and atlas...'); tic;
regopts = loadRegOpts(inputpath);
assert(isfield(regopts, 'affine_atlas_to_samp'), ...
    'multiobjCordRegistration:notInitialized', ...
    'Run initializeCordRegistration before multiobjCordRegistration.');

[tv, av, atlasinfo] = loadCordAtlasVolumes(regopts);
straightvol = readDownStack(regopts.straightvolpath);

refatlas    = imref3d(atlasinfo.regsize);
refsample   = imref3d(size(straightvol));
tvtemp      = medfilt3(tv);
usepointaff = false;
fprintf('Done! Took %2.2f s\n', toc);
%==========================================================================
% user control points
cppath        = dir(fullfile(regopts.savepath, 'corresponding_points.mat'));
cptshistology = zeros(0, 3);
cpaffine      = zeros(0, 3);
wtforpoints   = 0;
transaff      = regopts.affine_atlas_to_samp;

if ~isempty(cppath) && control_point_wt > 0
    cpdata = load(fullfile(cppath.folder, cppath.name));

    % make sure control points work
    np1           = cellfun(@(x) size(x,1), cpdata.atlas_control_points);
    np2           = cellfun(@(x) size(x,1), cpdata.histology_control_points);
    ikeep         = np1 == np2 & np1 > 0;
    cptsatlas     = cat(1, cpdata.atlas_control_points{ikeep});
    cptshistology = cat(1, cpdata.histology_control_points{ikeep});

    % make sure x is in the correct place
    cptsatlas     = cptsatlas(:, [2 1 3]);
    cptshistology = cptshistology(:, [2 1 3]);
    cpaffine      = cptsatlas;

    fprintf('Found %d user-defined control points. \n', size(cptshistology, 1));
    %----------------------------------------------------------------------
    % with enough of them, the affine transform is better refitted from the
    % points than kept from the image-based initialisation
    if size(cptshistology, 1) > 16
        atlasori    = regopts.affine_atlas_to_samp.transformPointsInverse(cptsatlas);
        transaff    = fitAffineTrans3D(atlasori, cptshistology);
        cpaffine    = transaff.transformPointsForward(atlasori);
        usepointaff = true;
    end
    %----------------------------------------------------------------------
    wtforpoints = control_point_wt;
    %----------------------------------------------------------------------
end
%==========================================================================
% warp atlas for the affine step
tvaffine = imwarp(tvtemp, refatlas, transaff, 'OutputView', refsample);
avaffine = imwarp(av, refatlas, transaff, 'nearest', 'OutputView', refsample);

volmax  = single(quantile(straightvol, 0.999, 'all'));
volplot = uint8(255 * single(straightvol) / volmax);

if usepointaff
    cf = plotCordAnnotation(volplot, avaffine);
    print(cf, fullfile(regopts.savepath, 'registration_point_affine'), '-dpng');
    close(cf);
end
%==========================================================================
% b-spline registration
[~, ~, bspltformpath, pathbspl] = performCordBsplineRegistration(...
    tvaffine, straightvol, regopts.registres(1) * 1e-3 * [1 1 1], ...
    cpaffine, cptshistology, wtforpoints, regopts.savepath, ...
    params.bspline_spatial_scale);

avreg = transformAnnotationVolume(bspltformpath, avaffine, 0.02);
%==========================================================================
cf = plotCordAnnotation(volplot, avreg);
print(cf, fullfile(regopts.savepath, 'registration_bspline'), '-dpng');
close(cf);
%==========================================================================
% invert it, so points measured on the sample can be pushed into the atlas
outdir    = fullfile(regopts.savepath, 'elastix_inverse_temp');
invstats  = invertElastixTransformCP(pathbspl, outdir);
tformpath = fullfile(regopts.savepath, ...
    sprintf('bspline_samp_to_atlas_%dum.txt', regopts.registres));
elastix_paramStruct2txt(tformpath, invstats.TransformParameters{1});
rmdir(invstats.outputDir, 's');
%==========================================================================
% the full chain from raw sample pixels to atlas voxels
params_pts_to_atlas = struct();
params_pts_to_atlas.samplekind       = 'cord';
params_pts_to_atlas.ori_pxsize       = regopts.pxsize;
params_pts_to_atlas.ori_size         = [regopts.Ny regopts.Nx regopts.Nz];
params_pts_to_atlas.registrationres  = regopts.registres;
params_pts_to_atlas.regvolsize       = regopts.regvolsize;
params_pts_to_atlas.how_to_perm      = regopts.sampleperm;
params_pts_to_atlas.samp_ikeepy      = regopts.yrange;
params_pts_to_atlas.samp_ikeepx      = regopts.xrange;
params_pts_to_atlas.samp_ikeeplong   = regopts.ikeeprange;
params_pts_to_atlas.straighten       = regopts.straighten;
params_pts_to_atlas.slicetforms      = regopts.slicetforms;
params_pts_to_atlas.tofliprc         = regopts.tofliprc;
params_pts_to_atlas.atlasres         = atlasinfo.atlasres;
params_pts_to_atlas.atlassize        = atlasinfo.regsize;
params_pts_to_atlas.atlassizenative  = atlasinfo.nativesize;
params_pts_to_atlas.elastix_um_to_mm = 1e-3;

params_pts_to_atlas.tform_bspline_samp20um_to_atlas_20um_px = tformpath;
params_pts_to_atlas.tform_affine_samp20um_to_atlas_20um_px  = transaff.invert;
params_pts_to_atlas.contol_pt_weight = control_point_wt;
params_pts_to_atlas.bspline_spatial_scale = params.bspline_spatial_scale;

save(fullfile(regopts.savepath, 'transform_params.mat'), '-struct', 'params_pts_to_atlas')
%==========================================================================
end

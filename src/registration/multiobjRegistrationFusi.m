function params_pts_to_atlas = multiobjRegistrationFusi(opts, contol_point_wt, usemultistep)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
opts.mousename = getOr(opts, 'mousename', 'x');
%==========================================================================
% read reg volume and control points
fprintf('Loading data volume and pre-computed clouds...'); tic;
% 
% load volume and control points
permvec    = opts.permute_sample_to_atlas;
% we permute the volume to match atlas
volume     = permuteBrainVolume(opts.sample, permvec);
regvolsize = size(volume);

facresize = opts.sample_res(abs(permvec))./opts.atlas_res;
volume    = imresize3(volume, 'Scale', facresize);

fprintf('Done! Took %2.2f s\n', toc);
%==========================================================================
optsreg.usemultistep              = usemultistep;
optsreg.n_histogram_bins          = getOr(opts, 'n_histogram_bins', 48);
optsreg.cpwt                      = contol_point_wt;
optsreg.bspline_spatial_scale     = getOr(opts, 'bspline_spatial_scale', 2);
opts.registres                    = opts.atlas_res;
optsreg.custom_sampleregion       = true;
%==========================================================================
% we also load user points
cppath        = dir(fullfile(opts.savepath,'*minimal.mat'));
cptshistology = zeros(0, 3);
cptsatlas     = zeros(0, 3);
if ~isempty(cppath)
    dpcp     = fullfile(cppath.folder,   cppath.name);
    cpdata   = load(dpcp);
    
    % make sure control points work
    np1           = cellfun(@(x) size(x,1),cpdata.atlas_points);
    np2           = cellfun(@(x) size(x,1),cpdata.sample_points);
    ikeep         = np1==np2 & np1>0;
    cptsatlas     = cat(1, cpdata.atlas_points{ikeep});
    cptshistology = cat(1, cpdata.sample_points{ikeep});
    cptsatlas     = cptsatlas(:, 1:3);
    cptshistology = cptshistology(:, 1:3);
    cptshistology = cptshistology.*facresize; % resize for isotropic

    % make sure x is in the correct place
    cptsatlas     = cptsatlas(:, [2 1 3]);
    cptshistology = cptshistology(:, [2 1 3]);
    fprintf('Found %d user-defined control points.\n', size(cptshistology, 1));    
end
%==========================================================================
% load atlas
tv   = opts.atlas;
av   = opts.annotation;
av   = resampleAnnotation(av, opts.parcelinfo, 'structure');
tv   = tv - imgaussfilt3(tv, 4*round( 0.3./opts.atlas_res )+1);
tv(tv<0) = 0;
tv = tv./quantile(tv,0.99,'all');
%==========================================================================
% affine optimization
% tform_aff       = fitAffineTrans3D(cptsatlas, cptshistology);
elastixparams = struct('NumberOfResolutions', 2, 'ImagePyramidSchedule', [2 1], ...
    'MaximumNumberOfIterations', [1000 1000],'NumberOfHistogramBins',optsreg.n_histogram_bins );
[~,~,affelastix] = performMultObjAffineRegistration(tv, volume, 1, ...
    cptsatlas, cptshistology, contol_point_wt*0.5, opts.savepath, elastixparams);
tform_aff = parse_elastix_tform_all(affelastix);
cpaffine        = tform_aff.transformPointsForward(cptsatlas);
%==========================================================================
fprintf('Warping Allen atlas... \n'); tic;
Rmoving  = imref3d(size(tv));
Rfixed   = imref3d(size(volume));
tvaffine = imwarp(tv, Rmoving, tform_aff, 'OutputView',Rfixed);
avaffine = imwarp(av, Rmoving, tform_aff, 'nearest','OutputView',Rfixed);
fprintf('Done! Took %2.2f s\n', toc);
%==========================================================================
% we plot the affine step
hilow = single(quantile(volume(volume>0), [0.05 0.96], 'all'));
voltoshow = uint8(255*(single(volume) - hilow(1))/range(hilow));
for idim = 1:3
    cf = plotAnnotationComparison(voltoshow, single(avaffine), idim, 'MarkerSize', 2);
    print(cf, fullfile(opts.savepath, sprintf('%s_dim%d_affine_registration', opts.mousename, idim)), '-dpng');
    close(cf);
end
%==========================================================================
%%
% we perform the b-spline registration
[reg, ~, bspltformpath, pathbspl] = fusiMultObjBsplineRegistration(tvaffine, volume, opts.registres, ...
    cpaffine, cptshistology, opts.savepath, optsreg);

avreg = transformAnnotationVolume(bspltformpath, avaffine, opts.registres);
%==========================================================================
% we plot the b-spline step
for idim = 1:3
    cf = plotAnnotationComparison(voltoshow, single(avreg), idim, 'MarkerSize', 2);
    print(cf, fullfile(opts.savepath, sprintf('%s_dim%d_bspline_registration', opts.mousename, idim)), '-dpng');
    close(cf);
end
%%
%%
%==========================================================================
% here we obtain the inverse transform
outdir     = fullfile(opts.savepath, 'elastix_inverse_temp');
invstats   = invertElastixTransformCP( pathbspl, outdir, 1, 1, 1);
tformpath  = fullfile(opts.savepath, 'bspline_samp_to_atlas_20um.txt');
elastix_paramStruct2txt(tformpath, invstats.TransformParameters{1});
rmdir(invstats.outputDir, 's'); % remove inversion directory
%==========================================================================
% data saving
params_pts_to_atlas = struct();
params_pts_to_atlas.atlasres         = opts.atlas_res;
params_pts_to_atlas.regvolsize       = regvolsize;
params_pts_to_atlas.atlassize        = size(tv);
params_pts_to_atlas.ori_pxsize       = opts.sample_res;
params_pts_to_atlas.ori_size         = size(opts.sample);
params_pts_to_atlas.how_to_perm      = opts.permute_sample_to_atlas;
params_pts_to_atlas.elastix_um_to_mm = 1e-3;

params_pts_to_atlas.tform_bspline_samp20um_to_atlas_20um_px = tformpath;
params_pts_to_atlas.tform_affine_samp20um_to_atlas_10um_px  = tform_aff.invert;
params_pts_to_atlas.contol_pt_weight = contol_point_wt;
params_pts_to_atlas.use_multistep    = usemultistep;
save(fullfile(opts.savepath, 'transform_params.mat'), '-struct', 'params_pts_to_atlas')
%==========================================================================
end
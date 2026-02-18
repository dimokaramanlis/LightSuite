function params_pts_to_atlas = multiobjCordRegistration(opts, contol_point_wt)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%==========================================================================
% read reg volume and control points
fprintf('Loading data volume and pre-computed clouds...'); tic;
optspath = dir(fullfile(opts.lsfolder, '*regopts.mat'));

% load volume and control points
dpopts    = fullfile(optspath.folder,   optspath.name);
regopts   = load(dpopts);
refatlas  = imref3d(size(regopts.tv));

straightvol = regopts.straightvol;
refsample   = imref3d(size(straightvol));
tvtemp      = medfilt3(regopts.tv);

% we permute the volume to match atlas
fprintf('Done! Took %2.2f s\n', toc);
%==========================================================================
% we also load user points
cppath        = dir(fullfile(opts.lsfolder,'corresponding_points.mat'));
cptshistology = zeros(0, 3);
cpaffine      = zeros(0, 3);
wtforpoints   = 0;
transaff      = regopts.affine_atlas_to_samp;

if ~isempty(cppath)
    dpcp     = fullfile(cppath.folder,   cppath.name);
    cpdata   = load(dpcp);
    
    % make sure control points work
    np1           = cellfun(@(x) size(x,1),cpdata.atlas_control_points);
    np2           = cellfun(@(x) size(x,1),cpdata.histology_control_points);
    ikeep         = np1==np2 & np1>0;
    cptsatlas     = cat(1, cpdata.atlas_control_points{ikeep});
    cptshistology = cat(1, cpdata.histology_control_points{ikeep});
    
    % make sure x is in the correct place
    cptsatlas     = cptsatlas(:, [2 1 3]);
    cptshistology = cptshistology(:, [2 1 3]);

    fprintf('Found %d user-defined control points. The affine transform will be refined.\n', size(cptshistology, 1)); 
    %-----------------------------------------------------------------------
    % we first refine the affine transform if the user has added points
    
    atlasuse = imwarp(tvtemp, refatlas, regopts.affine_atlas_to_samp, 'OutputView', refsample);
    
    [~,~, tformpath, ~]  = ...
        performMultObjAffineRegistration(atlasuse, straightvol, 1,...
        cptsatlas, cptshistology, contol_point_wt, regopts.lsfolder);
    newtrans = affinetform3d(parse_elastix_tform(tformpath));
    % transaff = fitAffineTrans3D(cptsatlas, cptshistology);
    
    transaff = affinetform3d(regopts.affine_atlas_to_samp.A*newtrans.A);
    cpaffine = newtrans.transformPointsForward(cptsatlas);
    %-----------------------------------------------------------------------
    wtforpoints = contol_point_wt;
    %-----------------------------------------------------------------------
end
%==========================================================================
% warp atlas for the affine step

tvaffine = imwarp(tvtemp, refatlas, transaff, 'OutputView',refsample);
avaffine = imwarp(regopts.av, refatlas, transaff, 'nearest', 'OutputView', refsample);

% we plot the affine step
volmax  = single(quantile(straightvol,0.999,'all'));
volplot = uint8(255*single(straightvol)/volmax);

cf     = plotCordAnnotation(volplot, avaffine);
print(cf, fullfile( regopts.lsfolder, 'registration_point_affine'), '-dpng')
close(cf);
%==========================================================================
% we perform the b-spline registration

[reg, ~, bspltformpath, pathbspl] = performCordBsplineRegistration(...
    tvaffine, straightvol, regopts.sampleres*1e-3, ...
    cpaffine, cptshistology, wtforpoints, regopts.lsfolder);

avreg = transformAnnotationVolume(bspltformpath, avaffine, 0.02);
%==========================================================================
cf     = plotCordAnnotation(volplot, avreg);
print(cf, fullfile( regopts.lsfolder, 'registration_bspline'), '-dpng')
close(cf);
%==========================================================================
% here we obtain the inverse transform
outdir     = fullfile(regopts.lsfolder, 'elastix_inverse_temp');
invstats   = invertElastixTransformCP( pathbspl, outdir);
tformpath  = fullfile(regopts.lsfolder, 'bspline_samp_to_atlas_20um.txt');
elastix_paramStruct2txt(tformpath, invstats.TransformParameters{1});
%==========================================================================
% data saving

params_pts_to_atlas = struct();
params_pts_to_atlas.tform_bspline_samp20um_to_atlas_20um_px = tformpath;
params_pts_to_atlas.tform_affine_samp20um_to_atlas_20um_px  = transaff.invert;
params_pts_to_atlas.contol_pt_weight                        = contol_point_wt;
params_pts_to_atlas.samp_ikeeplong = regopts.ikeeprange;
params_pts_to_atlas.samp_ikeepx    = regopts.xrange;
params_pts_to_atlas.samp_ikeepy    = regopts.yrange;
params_pts_to_atlas.how_to_perm    = regopts.sampleperm;
params_pts_to_atlas.slicetforms    = regopts.slicetforms;
params_pts_to_atlas.sampleres      = regopts.sampleres;
params_pts_to_atlas.tofliprc       = regopts.tofliprc;
params_pts_to_atlas.atlassize      = size(regopts.tv);

save(fullfile(regopts.lsfolder, 'transform_params.mat'), '-struct', 'params_pts_to_atlas')
%==========================================================================
end
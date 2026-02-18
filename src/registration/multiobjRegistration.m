function params_pts_to_atlas = multiobjRegistration(opts, contol_point_wt, usemultistep)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%==========================================================================
% read reg volume and control points
fprintf('Loading data volume and pre-computed clouds...'); tic;
optspath = dir(fullfile(opts.savepath,'regopts.mat'));

% load volume and control points
dp       = opts.regvolpath;
dpopts   = fullfile(optspath.folder,   optspath.name);

volume   = readDownStack(dp);
regopts  = load(dpopts);
regopts  = regopts.opts;

% we permute the volume to match atlas
volume  = permute(volume, regopts.permute_sample_to_atlas);
fprintf('Done! Took %2.2f s\n', toc);
%==========================================================================
% we load automated control points
autocptshistology = regopts.autocpsample;
autocptsatlas     = regopts.autocpatlas/regopts.downfac_reg;
distancethin      = 1000/regopts.atlasres;
if ~opts.augmentpoints
    distancethin = 3*distancethin;
end
%==========================================================================
% we also load user points
cppath        = dir(fullfile(opts.savepath,'*tform.mat'));
cptshistology = zeros(0, 3);
cptsatlas     = zeros(0, 3);
if ~isempty(cppath)
    dpcp     = fullfile(cppath.folder,   cppath.name);
    cpdata   = load(dpcp);
    
    % make sure control points work
    np1           = cellfun(@(x) size(x,1),cpdata.atlas_control_points);
    np2           = cellfun(@(x) size(x,1),cpdata.histology_control_points);
    ikeep         = np1==np2 & np1>0;
    cptsatlas     = cat(1, cpdata.atlas_control_points{ikeep});
    cptshistology = cat(1, cpdata.histology_control_points{ikeep});
    cptsatlas     = cptsatlas(:, 1:3);
    cptshistology = cptshistology(:, 1:3);


    % make sure x is in the correct place
    cptsatlas     = cptsatlas(:, [2 1 3]);
    cptshistology = cptshistology(:, [2 1 3]);
    cptshistology = regopts.original_trans.transformPointsInverse(cptshistology);
    cptsatlas     = cptsatlas/regopts.downfac_reg;

    Dall              = pdist2(cptsatlas, autocptsatlas);
    ikeepatlas        = all(Dall > distancethin, 1);

    autocptshistology = autocptshistology(ikeepatlas, :);
    autocptsatlas     = autocptsatlas(ikeepatlas, :);
    fprintf('Found %d user-defined control points.\n', size(cptshistology, 1));    
end
%==========================================================================
idxkeep         = thinPointList(autocptsatlas, distancethin);
fprintf('Augmentation with %d auto-extracted control points.\n', nnz(idxkeep));
afcptsatlas     = [cptsatlas; autocptsatlas(idxkeep, :)];
afcptshistology = [cptshistology; autocptshistology(idxkeep, :)];
tform_aff       = fitAffineTrans3D(afcptsatlas, afcptshistology);
cpaffine        = tform_aff.transformPointsForward(cptsatlas);
if size(cpaffine, 1) == 0
    % revert to automated mode
    contol_point_wt = contol_point_wt/2;
    cptshistology   = autocptshistology(idxkeep, :);
    cpaffine        = tform_aff.transformPointsForward(autocptsatlas(idxkeep, :));
end
%==========================================================================
% subplot(1,2,1)
% plot(cptshistology(:),cptsatlas(:),'o')
% title('Conrol points -  no alignment')
% subplot(1,2,2)
% plot(cptshistology(:),cpaffine(:),'o')
% title('Conrol points - after affine transform')

%==========================================================================
% load atlas
fprintf('Loading and warping Allen atlas... \n'); tic;
allen_atlas_path = fileparts(which('average_template_10.nii.gz'));
tv      = niftiread(fullfile(allen_atlas_path,'average_template_10.nii.gz'));
av      = niftiread(fullfile(allen_atlas_path,'annotation_10.nii.gz'));

Rmoving  = imref3d(size(tv));
Rfixed   = imref3d(size(volume));
tvaffine = imwarp(tv, Rmoving, tform_aff, 'OutputView',Rfixed);
avaffine = imwarp(av, Rmoving, tform_aff, 'nearest','OutputView',Rfixed);
fprintf('Done! Took %2.2f s\n', toc);
%==========================================================================
% we plot the affine step
voltoshow = uint8(255*single(volume)/single(quantile(volume, 0.999, 'all')));
for idim = 1:3
    cf = plotAnnotationComparison(voltoshow, single(avaffine), idim);
    print(cf, fullfile(regopts.savepath, sprintf('%s_dim%d_affine_registration', opts.mousename, idim)), '-dpng');
    close(cf);
end
%==========================================================================
%%
% we perform the b-spline registration

% tvaffine = single(tvaffine);
% tvaffine = tvaffine/quantile(tvaffine,0.999,'all');
% tvaffine = uint8(255 * tvaffine);
% 
% volume = single(volume);
% volume = volume/quantile(volume,0.999,'all');
% volume = uint16((2^16-1) * volume);%volume = uint8(255 * volume); %
%%
[reg, ~, bspltformpath, pathbspl] = performMultObjBsplineRegistration(tvaffine, volume, opts.registres*1e-3, ...
    cpaffine, cptshistology, contol_point_wt, usemultistep, regopts.savepath);

avreg = transformAnnotationVolume(bspltformpath, avaffine, opts.registres*1e-3);
%==========================================================================

% we plot the b-spline step
for idim = 1:3
    cf = plotAnnotationComparison(voltoshow, single(avreg), idim);
    print(cf, fullfile(regopts.savepath, sprintf('%s_dim%d_bspline_registration', opts.mousename, idim)), '-dpng');
    close(cf);
end
%%
%==========================================================================
% here we obtain the inverse transform
outdir     = fullfile(regopts.savepath, 'elastix_inverse_temp');
invstats   = invertElastixTransformCP( pathbspl, outdir);
tformpath  = fullfile(regopts.savepath, 'bspline_samp_to_atlas_20um.txt');
elastix_paramStruct2txt(tformpath, invstats.TransformParameters{1});
rmdir(invstats.outputDir, 's'); % remove inversion directory
%==========================================================================
% data saving

params_pts_to_atlas = struct();
params_pts_to_atlas.atlasres         = regopts.atlasres;
params_pts_to_atlas.regvolsize       = regopts.regvolsize;
params_pts_to_atlas.atlassize        = size(tv);
params_pts_to_atlas.ori_pxsize       = regopts.pxsize;
params_pts_to_atlas.ori_size         = [regopts.Ny regopts.Nx regopts.Nz];
params_pts_to_atlas.how_to_perm      = regopts.permute_sample_to_atlas;
params_pts_to_atlas.elastix_um_to_mm = 1e-3;

params_pts_to_atlas.tform_bspline_samp20um_to_atlas_20um_px = tformpath;
params_pts_to_atlas.tform_affine_samp20um_to_atlas_10um_px  = tform_aff.invert;
params_pts_to_atlas.contol_pt_weight = contol_point_wt;
params_pts_to_atlas.use_multistep    = usemultistep;
save(fullfile(opts.savepath, 'transform_params.mat'), '-struct', 'params_pts_to_atlas')
%==========================================================================
end


% viewerUnregistered = viewer3d(BackgroundColor="black",BackgroundGradient="off");
% volshow(volume,Parent=viewerUnregistered,RenderingStyle="Isosurface", ...
%     Colormap=[1 0 1],Alphamap=1);
% volshow(tvaffine,Parent=viewerUnregistered,RenderingStyle="Isosurface", ...
%     Colormap=[0 1 0],Alphamap=1);
% %
%%

%==========================================================================
% transfile = dir(fullfile(pathout, 'TransformParameters*.txt'));
% regpoints = transformix(cpaffine*0.02,fullfile(transfile.folder, transfile.name));
% msecurr   = mean(sqrt(sum((regpoints.OutputPoint - cptshistology*0.02).^2,2)));
% mseaffine = mean(sqrt(sum((cpaffine*0.02 - cptshistology*0.02).^2,2)));
% 
% fprintf('Control point distance is %d um\n', round(msecurr*1e3))
% voltoshow = uint8(255*single(volume)/180);
% plotSliceComparison(voltoshow, reg, 1)
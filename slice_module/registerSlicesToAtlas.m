function regparams = registerSlicesToAtlas(opts)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

regparams = struct();
%==========================================================================
% read reg volume and control points
volpath = dir(fullfile(opts.procpath,'*20um.tif'));
cppath   = dir(fullfile(opts.procpath,'*tform.mat'));
optspath = dir(fullfile(opts.procpath,'*regopts.mat'));

% load volume and control points
dp       = fullfile(volpath.folder, volpath.name);
dpcp     = fullfile(cppath.folder,   cppath.name);
dpopts   = fullfile(optspath.folder,   optspath.name);

% volume   = readDownStack(fullfile(opts.procpath, 'volume_for_inspection.tiff'), 1);

volume   = readDownStack(dp);
cpdata   = load(dpcp);
regopts  = load(dpopts);

% we permute the volume to match atlas
volume  = permute(volume, regopts.howtoperm);
Nslices = size(volume, 1);
%==========================================================================
downfac_reg = regopts.allenres/regopts.registres;

allen_atlas_path = fileparts(which('average_template_10.nii.gz'));
tv      = niftiread(fullfile(allen_atlas_path,'average_template_10.nii.gz'));
av      = niftiread(fullfile(allen_atlas_path,'annotation_10.nii.gz'));
tvdown  = imresize3(tv, downfac_reg);
avdown  = imresize3(av, downfac_reg, "Method", "nearest");


Ratlas         = imref3d(size(tvdown));
Rvolume        = imref3d(size(volume), 1, regopts.pxsizes(1), 1);
[tvnew, rnew]  = imwarp(tvdown, Ratlas, opts.tform_rigid_AllenToSample_20um, 'linear');
[avnew, rnew]  = imwarp(avdown, Ratlas, opts.tform_rigid_AllenToSample_20um, 'nearest');


xworld     = [Rvolume.XWorldLimits(1) + 0.5, Rvolume.XWorldLimits(2) - 0.5];
zworld     = [Rvolume.ZWorldLimits(1) + 0.5, Rvolume.ZWorldLimits(2) - 0.5];
zworld(1)  = max(zworld(1), rnew.ZWorldLimits(1));
zworld(2)  = min(zworld(2), rnew.ZWorldLimits(2));
yworld     = [Rvolume.YWorldLimits(1) + 0.5 - 10*opts.pxsizes(1),...
    Rvolume.YWorldLimits(2) + 10*opts.pxsizes(1) - 0.5];
yworld(1)  = max(yworld(1), rnew.YWorldLimits(1));
yworld(2)  = min(yworld(2), rnew.YWorldLimits(2));

[yy,xx,zz]     = rnew.worldToSubscript(xworld, yworld, zworld);
tvnewred       = tvnew(yy(1):yy(2), xx(1):xx(2), zz(1):zz(2));
avnew          = avnew(yy(1):yy(2), xx(1):xx(2), zz(1):zz(2));


% yatlasvals     = linspace(yworld(1), yworld(2), yy(2)-yy(1)+2);
% yatlasvals     = yatlasvals(1:end-1) + median(diff(yatlasvals))/2;
% ysamplevals    = linspace(Rvolume.YWorldLimits(1), Rvolume.YWorldLimits(2), Nslices+1);
% ysamplevals    = ysamplevals(1:end-1) + median(diff(ysamplevals))/2;
% [~, atlasinds] = min(pdist2(ysamplevals',yatlasvals'), [],2);


atlasinds       = cellfun(@(x) x(1,1), cpdata.atlas_control_points);

%%
% clf;
% scatter3(cptshistology(:, 1), cptshistology(:, 2), cptshistology(:, 3))
% hold on;
% scatter3(cpatlasori(:, 1), cpatlasori(:, 2), cpatlasori(:, 3))
% 
% tformout = fitRigidTrans3D(cpatlasori, cptshistology);
% % tformout = fitAffineTrans3D(cpatlasori, cptshistology);
% 
% newatlas = tformout.transformPointsForward(cpatlasori);
% clf;
% scatter3(cptshistology(:, 1), cptshistology(:, 2), cptshistology(:, 3))
% hold on;
% scatter3(newatlas(:, 1), newatlas(:, 2), newatlas(:, 3))
% 
% raatlas  = imref3d(size(tv),  0.5, 0.5, 0.5);
% rasample = imref3d(size(volume), 1, regopts.pxsizes(1), 1);
% tvtest   = imwarp(tv, raatlas, tformout,'linear', 'OutputView',rasample);
% tvoritest = imwarp(tv, raatlas, opts.tform_rigid_AllenToSample_20um,'linear', 'OutputView',rasample);
%%
% clf;
% for ishow = 1:Nslices
% % ishow = 45;
%     subplot(1,3,1)
%     imagesc(squeeze(volume(ishow, :,:))); axis equal; axis tight;
%     subplot(1,3,2)
%     imagesc(squeeze(tvoritest(ishow,:,:))); axis equal; axis tight;
%     subplot(1,3,3)
%     imagesc(squeeze(tvnew(atlasinds(ishow),:,:))); axis equal; axis tight;
%     pause;
% end
%%
%=========================================================================
% for every slice, we fit an affine transform from atlas to the slice
% if no control points are available, we use point clouds
fprintf('Using control points and elastix atlas fitting...\n'); 
wholetic = tic;
msg = repmat('=', [1 100]);
ratlas      = imref2d(size(tvnewred,  [2 3]));
rahist      = imref2d(size(volume, [2 3]));
slicetforms = affinetform2d;
cpwt        = 0.2;
forsavepath = fullfile(opts.procpath, 'elastix_forward');
revsavepath = fullfile(opts.procpath, 'elastix_reverse');
makeNewDir(forsavepath);
makeNewDir(revsavepath);

for islice = 1:Nslices
    %------------------------------------------------------------------
    fprintf('Aligning slice %d/%d to atlas using elastix...\n', islice, Nslices); 
    slicetimer   = tic;

    histim    = squeeze(volume(islice, :,:));
    atlasim   = squeeze(tvnewred(atlasinds(islice),:,:));
    annotim   = squeeze(avnew(atlasinds(islice),:,:));
    
    fixedpts  = cpdata.histology_control_points{islice}(:, [3 2]);
    movingpts = cpdata.atlas_control_points{islice}(:, [3 2]);
    %------------------------------------------------------------------
    % affine estimation

    if isempty(movingpts) || (size(movingpts, 1) < 5)
        % TO DO!!!!!
        tformcurr = fitPointCloudsAffine(atlasim, histim, opts.registres);
    else
        tformcurr = fitgeotform2d(movingpts, fixedpts, 'affine');
    end

    slicetforms(islice, 1) = tformcurr;
    %------------------------------------------------------------------
    % bspline estimation
    movptsaffine = tformcurr.transformPointsForward(movingpts);
    affatlasim   = imwarp(atlasim, ratlas, tformcurr, "linear",  "OutputView", rahist,...
        'FillValue', 1);
    affannotim = imwarp(annotim, ratlas, tformcurr, "nearest", "OutputView", rahist);
    dpsavefor  = fullfile(forsavepath, sprintf('%03d_slice', islice));
    [regimg, tformpath] = bsplineRegisterSlice(affatlasim, histim, opts.registres*1e-3, ...
        movptsaffine, fixedpts, cpwt, dpsavefor);

    % we rename the transform
    slicename     = sprintf('%03d_slice_bspline_atlas_to_samp_20um.txt', islice);
    bspltformpath = fullfile(forsavepath, slicename);
    copyfile(tformpath.TransformParametersFname{1}, bspltformpath);
    %------------------------------------------------------------------
    % transformix for illustration
    avreg     = transformAnnotationVolume(bspltformpath, affannotim, opts.registres*1e-3);
    sliceplot = uint8(255 * single(histim)/single(quantile(histim, 0.999, 'all')));
    cf = plotRegistrationComparison(sliceplot, cat(3, affannotim, avreg), ...
        {'affine (control points)', 'bspline'});
    savepngFast(cf, forsavepath, sprintf('%03d_slice_registration_comparison', islice), 300, 2);
    close(cf);
    %------------------------------------------------------------------
    % inverting
    dpsaverev    = fullfile(revsavepath, sprintf('%03d_slice', islice));
    invstats     = invertElastixTransformCP(dpsavefor, dpsaverev);
    slicenamerev = sprintf('%03d_slice_bspline_samp_to_atlas_20um.txt', islice);
    revtformpath = fullfile(revsavepath, slicenamerev);
    elastix_paramStruct2txt(revtformpath, invstats.TransformParameters{1});
    %------------------------------------------------------------------
    rmdir(dpsavefor, 's');
    rmdir(dpsaverev, 's');
    fprintf('Finished with slice! Took %2.2f s.\n%s\n', toc(slicetimer), msg); 
    %------------------------------------------------------------------
end
% we save the affine registrations to the output structure
regparams.affine_tform_atlas_to_image = slicetforms;
fprintf('DONE with all slices! Took %d min\n', ceil(toc(wholetic)/60));
%%
%=========================================================================
%%
% testvol = zeros([ratlas.ImageSize Nslices], 'uint16');
% for islice = 1:Nslices
%     histim  = squeeze(volume(islice, :,:));
%     atlasim = squeeze(tvnew(atlasinds(islice),:,:));
%     revtformpath = dir(fullfile(revsavepath, sprintf('%03d_slice_*', islice)));
%     revtformpath = fullfile(revtformpath.folder, revtformpath.name);
%     testim  = transformix(histim,revtformpath, 'movingscale', 0.01*[1 1], 'verbose', 0);
%     testim  = uint16(abs(testim));
%     testvol(:,:,islice)  = imwarp(testim, rahist, slicetforms(islice).invert, 'OutputView',ratlas);
% end
% % imshowpair(testim, atlasim)
% 
% 
% % we have to take the transformed atlas and zero it out
% % we fill in the atlas the slices we chose to register
% % we undo the rigid transform to go back to the original atlas space
% 
% rcurr = imref3d([1 size(testim)], 1,  7.5,1);
% [avnew, rnew]  = imwarp(reshape(testim, [1, size(testim)]),...
%     rcurr, opts.tform_rigid_AllenToSample_20um.invert, 'nearest');

%%
%=========================================================================
%%

%==========================================================================
% data saving

% params_pts_to_atlas = struct();
% params_pts_to_atlas.atlasres         = regopts.atlasres;
% params_pts_to_atlas.regvolsize       = regopts.regvolsize;
% params_pts_to_atlas.atlassize        = size(tv);
% params_pts_to_atlas.ori_pxsize       = regopts.pxsize;
% params_pts_to_atlas.ori_size         = [regopts.Ny regopts.Nx regopts.Nz];
% params_pts_to_atlas.how_to_perm      = regopts.permute_sample_to_atlas;
% params_pts_to_atlas.elastix_um_to_mm = 1e-3;
% 
% params_pts_to_atlas.tform_bspline_samp20um_to_atlas_20um_px = tformpath;
% params_pts_to_atlas.tform_affine_samp20um_to_atlas_10um_px  = tform_aff.invert;
% params_pts_to_atlas.contol_pt_weight = contol_point_wt;
% params_pts_to_atlas.use_multistep    = usemultistep;
save(fullfile(opts.procpath, 'transform_params.mat'), '-struct', 'regparams')
%==========================================================================
end
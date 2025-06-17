function regparams = registerSlicesToAtlas(opts)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

regparams = struct();
%==========================================================================
% read reg volume and control points
volpath  = dir(fullfile(opts.procpath,'*20um.tif'));
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
nfac        = regopts.extentfactor;


allen_atlas_path = fileparts(which('average_template_10.nii.gz'));
tv      = niftiread(fullfile(allen_atlas_path,'average_template_10.nii.gz'));
av      = niftiread(fullfile(allen_atlas_path,'annotation_10.nii.gz'));
tv      = tv(regopts.atlasaplims(1):regopts.atlasaplims(2), :, :);
av      = av(regopts.atlasaplims(1):regopts.atlasaplims(2), :, :);
tvdown  = imresize3(tv, downfac_reg);
avdown  = imresize3(av, downfac_reg, "Method", "nearest");

Ratlas  = imref3d(size(tvdown));
Rvolume = imref3d(size(volume), 1, regopts.pxsizes(1), 1);
yworld  = [Rvolume.YWorldLimits(1)-regopts.pxsizes(1)*nfac, Rvolume.YWorldLimits(2)+nfac*regopts.pxsizes(1)];
ypix    = range(yworld);
Rout    = imref3d([ypix, size(tvdown, [2 3])], Rvolume.XWorldLimits, yworld,Rvolume.ZWorldLimits);

[tvnew, rnew]  = imwarp(tvdown, Ratlas, opts.tformrigid_allen_to_samp_20um, 'linear',  'OutputView', Rout);
[avnew, rnew]  = imwarp(avdown, Ratlas, opts.tformrigid_allen_to_samp_20um, 'nearest', 'OutputView', Rout);

% get original prediction indices
yatlasvals       = linspace(yworld(1), yworld(2), ypix + 1);
yatlasvals       = yatlasvals(1:end-1) + median(diff(yatlasvals))/2;
ysamplevals      = linspace(Rvolume.YWorldLimits(1), Rvolume.YWorldLimits(2), Nslices+1);
ysamplevals      = ysamplevals(1:end-1) + median(diff(ysamplevals))/2;
[~, atlasinds]   = min(pdist2(ysamplevals',yatlasvals'), [],2);

hascp            = ~cellfun(@isempty, cpdata.atlas_control_points);
useratlasinds    = cellfun(@(x) x(1,1), cpdata.atlas_control_points(hascp));
if nnz(hascp) > 4
    % refine remaining
    pfit      = polyfit(atlasinds(hascp), useratlasinds, 1);
    atlasinds = round(polyval(pfit, atlasinds));
end
atlasinds(hascp) = useratlasinds;

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
ratlas      = imref2d(size(tvnew,  [2 3]));
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
    atlasim   = squeeze(tvnew(atlasinds(islice),:,:));
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
fprintf('DONE with all slices! Took %d min\n', ceil(toc(wholetic)/60));
%%
%=========================================================================




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

regparams.atlasres         = regopts.allenres;
regparams.how_to_perm      = regopts.howtoperm;
regparams.atlasaplims      = regopts.atlasaplims;
regparams.tformrigid_allen_to_samp_20um       = regopts.tformrigid_allen_to_samp_20um;
regparams.tformbspline_samp20um_to_atlas_20um = revsavepath;
regparams.tformaffine_tform_atlas_to_image    = slicetforms;
regparams.space_sample_20um                   = rnew;
regparams.space_atlas_20um                    = ratlas;
regparams.space3d_sample_20um                 = Rout;
regparams.space3d_atlas_20um                  = Ratlas;
regparams.sliceids_in_sample_space            = atlasinds;
save(fullfile(opts.procpath, 'transform_params.mat'), '-struct', 'regparams')
%==========================================================================
end
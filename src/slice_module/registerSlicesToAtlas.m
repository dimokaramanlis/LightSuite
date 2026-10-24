function regparams = registerSlicesToAtlas(opts)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

regparams = struct();
%==========================================================================
% read reg volume and options
volpath  = dir(fullfile(opts.procpath,'*20um.tif'));
optspath = dir(fullfile(opts.procpath,'*regopts.mat'));

dp       = fullfile(volpath.folder, volpath.name);
dpopts   = fullfile(optspath.folder,   optspath.name);

volume   = readDownStack(dp);
regopts  = load(dpopts);

% we permute the volume to match atlas
volume  = permute(volume, regopts.howtoperm);
Nslices = size(volume, 1);
%==========================================================================
% load control points if available
cppath   = dir(fullfile(opts.procpath,'*tform.mat'));
if ~isempty(cppath)
    dpcp              = fullfile(cppath.folder, cppath.name);
    cpdata            = load(dpcp);
    histology_cpoints = cpdata.histology_control_points;
    atlas_cpoints     = cpdata.atlas_control_points;
else
    histology_cpoints = repmat({zeros(0,2)}, [Nslices 1]);
    atlas_cpoints     = repmat({zeros(0,2)}, [Nslices 1]);
end
%==========================================================================
downfac_reg = regopts.allenres/regopts.registres;
nfac        = ceil(regopts.extentfactor * 7.5/regopts.pxsizes(1));


allen_atlas_path = fileparts(which('average_template_10.nii.gz'));
tv      = niftiread(fullfile(allen_atlas_path,'average_template_10.nii.gz'));
av      = niftiread(fullfile(allen_atlas_path,'annotation_10.nii.gz'));
tv      = tv(regopts.atlasaplims(1):regopts.atlasaplims(2), :, :);
av      = av(regopts.atlasaplims(1):regopts.atlasaplims(2), :, :);
tvdown  = imresize3(tv, downfac_reg);
avdown  = imresize3(av, downfac_reg, "Method", "nearest");

Ratlas  = imref3d(size(tvdown));
Rvolume = imref3d(size(volume), 1, regopts.pxsizes(1), 1);
%--------------------------------------------------------------------------
% Check if the cutting angle data file exists
tformuse        = regopts.tformrigid_allen_to_samp_20um;
angle_data_path = fullfile(opts.procpath, 'cutting_angle_data.mat');

if exist(angle_data_path, 'file')
    fprintf('Found saved cutting angle data. Overwriting transformation rotation...\n');
    % --- 1. Load and process angle data ---
    data = load(angle_data_path);
    cutting_angle_data = data.cutting_angle_data;
    tformuse = applyAngleToTransform(tformuse, Ratlas, cutting_angle_data);
end
% same output view as matchControlPointsInSlices (covers the whole atlas in AP)
[Rout, yworld, ypix, nlow] = getSliceAtlasOutputRef(Rvolume, Ratlas, tformuse, nfac, regopts.pxsizes(1));
%--------------------------------------------------------------------------

[tvnew, rnew]  = imwarp(tvdown, Ratlas, tformuse, 'linear',  'OutputView', Rout);
[avnew, rnew]  = imwarp(avdown, Ratlas, tformuse, 'nearest', 'OutputView', Rout);

% get original prediction indices
yatlasvals       = linspace(yworld(1), yworld(2), ypix + 1);
yatlasvals       = yatlasvals(1:end-1) + median(diff(yatlasvals))/2;
ysamplevals      = linspace(Rvolume.YWorldLimits(1), Rvolume.YWorldLimits(2), Nslices+1);
ysamplevals      = ysamplevals(1:end-1) + median(diff(ysamplevals))/2;
[~, atlasinds]   = min(pdist2(ysamplevals',yatlasvals'), [],2);

%=========================================================================
% deal with missing control points
hascp            = ~cellfun(@isempty,   atlas_cpoints);
useratlasinds    = cellfun(@(x) x(1,1), atlas_cpoints(hascp)) + nlow; % saved in legacy AP range
if nnz(hascp) > 3
    % refine remaining
    % pfit      = polyfit(atlasinds(hascp), useratlasinds, 1);
    % atlasinds = round(polyval(pfit, atlasinds));
    ivals     = interp1(find(hascp), useratlasinds, 1:Nslices, 'linear', 'extrap')';
    atlasinds = round(ivals);
end
atlasinds(hascp) = useratlasinds;
atlasinds        = min(max(round(atlasinds), 1), ypix);
%=========================================================================
% for every slice, we fit an affine transform from atlas to the slice
% if no control points are available, we use only image info
fprintf('Using control points and elastix atlas fitting...\n'); 

wholetic = tic; msg = repmat('=', [1 100]);

ratlas      = imref2d(size(tvnew,  [2 3]));
rahist      = imref2d(size(volume, [2 3]));
slicetforms = affinetform2d;
cpwt        = getOr(opts, 'cpwt',0.2);
% B-spline grid spacing in mm, how much you allow the bspline to bend
bsplscale   = getOr(opts, 'bspline_spatial_scale', 0.96);
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
    
    fixedpts  = histology_cpoints{islice}(:, [3 2]);
    fixedpts  = reshape(fixedpts, [], 2);
    movingpts = atlas_cpoints{islice}(:, [3 2]);
    movingpts = reshape(movingpts, [], 2);
    Nmov      = size(movingpts, 1);
    %------------------------------------------------------------------
    % timeshistology = histology_cpoints{islice}(:, 4);
    % timesatlas     = atlas_cpoints{islice}(:, 4);
    %------------------------------------------------------------------
    % affine estimation

    if isempty(movingpts) || (Nmov < 5)
        tformcurr = fitPointCloudsAffine(atlasim, histim, opts.registres);
        tstruse   = sprintf('Not enough control points, using only image...\n\n');
    else
        tformcurr = fitgeotform2d(movingpts, fixedpts, 'affine');
        tstruse   = sprintf('%d control points found, thanks for the effort!\n\n', Nmov);
    end

    slicetforms(islice, 1) = tformcurr;
    fprintf(tstruse)
    %------------------------------------------------------------------
    % bspline estimation
    movptsaffine = tformcurr.transformPointsForward(movingpts);
    affatlasim   = imwarp(atlasim, ratlas, tformcurr, "linear",  "OutputView", rahist,...
        'FillValue', 1);
    affannotim = imwarp(annotim, ratlas, tformcurr, "nearest", "OutputView", rahist);
    dpsavefor  = fullfile(forsavepath, sprintf('%03d_slice', islice));
    % grid spacing adapts to how well the points cover the brain in the slice
    [gridspacing, coarsescale, meddist] = getAdaptiveBsplineScale(affannotim > 0, ...
        fixedpts, opts.registres*1e-3, bsplscale);
    fprintf('B-spline grid spacing: %2.2f mm (coarse = %2.2f mm, median distance to points = %2.2f mm)\n', ...
        gridspacing, coarsescale, meddist);
    [regimg, tformpath] = bsplineRegisterSlice(affatlasim, histim, opts.registres*1e-3, ...
        movptsaffine, fixedpts, cpwt, dpsavefor, gridspacing);

    % we rename the transform
    slicename     = sprintf('%03d_slice_bspline_atlas_to_samp_20um.txt', islice);
    bspltformpath = fullfile(forsavepath, slicename);
    copyfile(tformpath.TransformParametersFname{1}, bspltformpath);
    %------------------------------------------------------------------
    % transformix for illustration
    avreg     = transformAnnotationVolume(bspltformpath, affannotim, opts.registres*1e-3);
    sliceplot = single(histim);
    slicelims  = quantile(sliceplot(sliceplot>0), [0.01 0.99], 'all');
    sliceplot = uint8(255 * (sliceplot - slicelims(1))/range(slicelims));
    txtstr1   = sprintf('affine (Npts = %d)', Nmov);
    cf = plotRegistrationComparison(sliceplot, cat(3, affannotim, avreg), ...
        {txtstr1, sprintf('bspline (scale = %2.2f mm)', gridspacing)}, fixedpts);
    print(cf, fullfile( forsavepath, sprintf('%03d_slice_registration_comparison', islice)), '-dpng')
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
%=========================================================================
% data saving

regparams.atlasres         = regopts.allenres;
regparams.how_to_perm      = regopts.howtoperm;
regparams.atlasaplims      = regopts.atlasaplims;
regparams.processres       = regopts.processres;
regparams.registres        = regopts.registres;
regparams.tformrigid_allen_to_samp_20um       = tformuse;
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

function [gridspacing, coarsescale, meddist] = getAdaptiveBsplineScale(brainmask, pts, pxsize, finescale)
% GETADAPTIVEBSPLINESCALE picks the B-spline grid spacing (mm) for one slice.
%   The coarsest spacing is half the sqrt of the brain area. With control
%   points, the spacing follows the median distance from brain pixels to the
%   nearest point, so clustered points barely refine the grid while
%   well-spread points bring it down to FINESCALE.
%   PTS are [x y] pixel coordinates in the same frame as BRAINMASK.
%==========================================================================
if ~any(brainmask, 'all')
    brainmask = true(size(brainmask));
end
coarsescale = 0.5 * sqrt(nnz(brainmask) * pxsize^2);

% keep points that fall inside the image
pts   = round(pts);
inimg = all(pts >= 1, 2) & pts(:,1) <= size(brainmask, 2) & pts(:,2) <= size(brainmask, 1);
pts   = pts(inimg, :);

if isempty(pts)
    meddist     = Inf;
    gridspacing = coarsescale;
    return;
end

ptmask = false(size(brainmask));
ptmask(sub2ind(size(brainmask), pts(:,2), pts(:,1))) = true;
distmm  = bwdist(ptmask) * pxsize;
meddist = median(distmm(brainmask));

gridspacing = min(max(meddist, finescale), coarsescale);
end

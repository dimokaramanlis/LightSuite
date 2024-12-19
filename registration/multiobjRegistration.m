function params_pts_to_atlas = multiobjRegistration(opts, contol_point_wt)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%==========================================================================
% read reg volume and control points
volpath = dir(fullfile(opts.savepath,'*.tif'));
cppath  = dir(fullfile(opts.savepath,'*tform.mat'));
optspath = dir(fullfile(opts.savepath,'*regopts.mat'));

% load volume and control points
dp      = fullfile(volpath.folder, volpath.name);
dpcp    = fullfile(cppath.folder,   cppath.name);
dpopts  = fullfile(optspath.folder,   optspath.name);

volume  = readDownStack(dp);
cpdata  = load(dpcp);
regopts = load(dpopts);
regopts = regopts.opts;

% make sure control points work
np1           = cellfun(@(x) size(x,1),cpdata.atlas_control_points);
np2           = cellfun(@(x) size(x,1),cpdata.histology_control_points);
ikeep         = np1==np2 & np1>0;
cptsatlas     = cat(1, cpdata.atlas_control_points{ikeep});
cptshistology = cat(1, cpdata.histology_control_points{ikeep});

% make sure x is in the correct place
cptsatlas     = cptsatlas(:, [2 1 3]);
cptshistology = cptshistology(:, [2 1 3]);
%==========================================================================
% 1st transform - AFFINE BASED ON CONTROL POINTS 
% TODO: Do affine with control points + intensity?
tform_aff     = fitAffineTrans3D(cptsatlas, cptshistology);
cpaffine      = tform_aff.transformPointsForward(cptsatlas);
% 
% subplot(1,2,1)
% plot(cptshistology(:),cptsatlas(:),'o')
% title('Conrol points -  no alignment')
% subplot(1,2,2)
% plot(cptshistology(:),cpaffine(:),'o')
% title('Conrol points - after affine transform')
%==========================================================================
% load atlas
fprintf('Loading and warping Allen atlas... \n'); tic;
allen_atlas_path = fileparts(which('template_volume_10um.npy'));
tv = readNPY(fullfile(allen_atlas_path,'template_volume_10um.npy'));
factv = 255/single(max(tv,[],"all"));
tv = uint8(single(tv)*factv);
tvdown = imresize3(tv,regopts.downfac_reg);
av = readNPY(fullfile(allen_atlas_path,'annotation_volume_10um_by_index.npy'));
av = imresize3(av,regopts.downfac_reg, "Method","nearest");


Rmoving  = imref3d(size(tvdown));
Rfixed   = imref3d(size(volume));
tvaffine = imwarp(tvdown, Rmoving, tform_aff, 'OutputView',Rfixed);
avaffine = imwarp(av, Rmoving, tform_aff, 'nearest','OutputView',Rfixed);
fprintf('Done! Took %2.2f s\n', toc);
%==========================================================================
% we plot the affine step
voltoshow = uint8(255*single(volume)/150);
for idim = 1:3
    cf = plotAnnotationComparison(voltoshow, single(avaffine), idim);
    savepngFast(cf, regopts.savepath, sprintf('dim%d_affine_registration', idim), 300, 2);
    close(cf);
end
%==========================================================================
%%
% we perform the b-spline registration
[reg, ~, bspltformpath, pathbspl] = performMultObjBsplineRegistration(tvaffine, volume, 0.02, ...
    cpaffine, cptshistology, contol_point_wt, regopts.savepath);
avreg = transformAnnotationVolume(bspltformpath, avaffine, 0.02);
%==========================================================================
%%
% we plot the b-spline step
for idim = 1:3
    cf = plotAnnotationComparison(voltoshow, single(avreg), idim);
    savepngFast(cf, regopts.savepath, sprintf('dim%d_bspline_registration', idim), 300, 2);
    close(cf);
end

%==========================================================================
% here we obtain the inverse transform
outdir     = fullfile(regopts.savepath, 'elastix_inverse_temp');
invstats   = invertElastixTransformCP( pathbspl, outdir);

tformpath = fullfile(regopts.savepath, 'bspline_samp_to_atlas_20um.txt');
copyfile(invstats.TransformParametersFname{1}, tformpath)
%==========================================================================
% data saving

params_pts_to_atlas = struct();
params_pts_to_atlas.atlasres         = regopts.atlasres;
params_pts_to_atlas.regvolsize       = regopts.regvolsize;
params_pts_to_atlas.regvolsizedwon   = regopts.regvolsize_down;
params_pts_to_atlas.atlassize        = size(tv);
params_pts_to_atlas.ori_pxsize       = regopts.pxsize;
params_pts_to_atlas.ori_size         = [regopts.Ny regopts.Nx regopts.Nz];
params_pts_to_atlas.how_to_perm      = regopts.permute_sample_to_atlas;
params_pts_to_atlas.elastix_um_to_mm = 1e-3;

params_pts_to_atlas.tform_rigid_samp20um_to_atlas_20um_px   = regopts.original_trans;
params_pts_to_atlas.tform_bspline_samp20um_to_atlas_20um_px = tformpath;
params_pts_to_atlas.tform_affine_samp20um_to_atlas_20um_px  = tform_aff.invert;

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
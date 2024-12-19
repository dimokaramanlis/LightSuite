
% atlas need to be set at the same dimension as the main thing
dp = 'D:\DATA_folder\Mice\DK031\Anatomy\vol_back.tif';
dpcp = 'D:\DATA_folder\Mice\DK031\Anatomy\atlas2histology_tform.mat';
volume    = readDownStack(dp);

volume    = permute(volume, [1 3 2]);

allen_atlas_path = fileparts(which('template_volume_10um.npy'));
tv = readNPY(fullfile(allen_atlas_path,'template_volume_10um.npy'));
av = readNPY(fullfile(allen_atlas_path,'annotation_volume_10um.npy'));
% [tv,av] = readGubraAtlas();

%%
viewerUnregistered = viewer3d(BackgroundColor="black",BackgroundGradient="off");
volshow(volume,Parent=viewerUnregistered,RenderingStyle="Isosurface", ...
    Colormap=[1 0 1],Alphamap=1);
volshow(tv,Parent=viewerUnregistered,RenderingStyle="Isosurface", ...
    Colormap=[0 1 0],Alphamap=1);
%%
% Options = struct();
cnew = tform.transformPointsForward(cptsatlas);
% 
% % Options.Spacing = [2^6 2^6 2^6];
% Options.Registration ='NonRigid';
% Options.Similarity   = 'mi';
% Options.Points1 = cnew;
% Options.Points2 = cptshistology;
% Options.PStrength=ones(size(cptsatlas,1),1);
% [Ireg,O_trans,Spacing,M,B,F] = image_registration(volwrap2, volume, Options);
% Opts.Verbose = 2;
% [Ireg,O_trans,Xreg] = point_registration(size(volume),cptshistology,cnew, Opts);
%%
melastixpath = 'C:\Users\karamanl\Documents\GitHub\matlab_elastix';
yamlmatlabpath = 'C:\Users\karamanl\Documents\GitHub\yamlmatlab';
addpath(genpath(melastixpath), genpath(yamlmatlabpath));

params = struct();
params.Transform = 'AffineTransform';
params.Optimizer = 'AdaptiveStochasticGradientDescent';
params.Metric    = 'AdvancedMattesMutualInformation';
params.MaximumNumberOfIterations=500;
params.NumberOfSpatialSamples=1E3;
params.FixedImageDimension  = 3;
params.MovingImageDimension = 3;
% params.BSplineInterpolationOrder = 1;
params.FinalBSplineInterpolationOrder = 3;
params.NumberOfResolutions  = 4;
params.SP_a     = 4000;
params.SP_A     = 20;
params.SP_alpha = 0.6;
params.NumberOfHistogramBins = 32;

reg=elastix(tv,volume,'D:\lightsheet\elout','elastix_default.yml','paramstruct',params);
%%
viewerUnregistered = viewer3d(BackgroundColor="black",BackgroundGradient="off");
volshow(volume,Parent=viewerUnregistered,RenderingStyle="Isosurface", ...
    Colormap=[1 0 1],Alphamap=1);
volshow(reg,Parent=viewerUnregistered,RenderingStyle="Isosurface", ...
    Colormap=[0 1 0],Alphamap=1);

viewerUnregistered = viewer3d(BackgroundColor="black",BackgroundGradient="off");
volshow(volume,Parent=viewerUnregistered,RenderingStyle="Isosurface", ...
    Colormap=[1 0 1],Alphamap=1);
volshow(reg,Parent=viewerUnregistered,RenderingStyle="Isosurface", ...
    Colormap=[0 1 0],Alphamap=1);
%%
paramFiles = {'D:\lightsheet/TransformParameters.1.txt'};
% newav = transformix(volwrap,'D:\lightsheet\elout',1);
regpoints=transformix(cnew,'D:\lightsheet\elout',1);
plot(cptshistology(:), cnew(:),'o', cptshistology(:), regpoints.OutputPoint(:),'o',...
    [0 1]*1000,[0 1]*1000)

%%

islice= 50;
curr_slice_warp    = squeeze(newav(islice,:,:));
av_warp_boundaries = round(conv2(curr_slice_warp,ones(3)./9,'same')) ~= curr_slice_warp;


%%

islice =80;
ppanel = panel();
ppanel.pack('h',2)
ppanel(1).select();
imshowpair(squeeze(volume(islice,:,:)), squeeze(reg(islice,:,:)))
title('Affine (control-point-based)')
ppanel(2).select();
imshowpair(squeeze(volume(islice,:,:)), squeeze(reg(islice,:,:)))
title('B-spline (elastix)')
%%

%%
subplot(1,3,1)
imagesc(squeeze(volume(islice,:,:)))
axis equal; axis tight;
colormap(gray)

subplot(1,3,2)
imagesc(squeeze(volwrap2(islice,:,:)))
axis equal; axis tight;
colormap(gray)

subplot(1,3,3)
imagesc(squeeze(reg(islice,:,:)));
axis equal; axis tight;
colormap(gray)






%%
clf; plotBrainGrid; hold on;plot3(cptsatlas(:,2),cptsatlas(:,3),cptsatlas(:,1),'ro')
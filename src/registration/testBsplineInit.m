function [regimg,tform_bspline, tformpath, pathtemp] = testBsplineInit(movingvol,fixedvol,...
    volscale, initialtransform, savepath)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%==========================================================================
addElastixRepoPaths;
params = struct();
%==========================================================================
% general parameters, probably won't touch
params.Registration                  = 'MultiResolutionRegistration';
params.Metric                        = 'AdvancedMattesMutualInformation';
params.Transform                       = 'RecursiveBSplineTransform';%'RecursiveBSplineTransform';
params.Optimizer                       = 'AdaptiveStochasticGradientDescent';
params.ImageSampler                    = 'RandomCoordinate';
params.AutomaticParameterEstimation    = true;
params.AutomaticScalesEstimation       = true;
params.BSplineInterpolationOrder       = 3;
params.FinalBSplineInterpolationOrder  = 3;
params.FixedImageDimension             = 3;
params.MovingImageDimension            = 3;
params.FixedImagePyramid               = 'FixedRecursiveImagePyramid';
params.MovingImagePyramid              = 'MovingRecursiveImagePyramid';
params.UseRandomSampleRegion           = true;
params.NewSamplesEveryIteration        = true;
params.NumberOfResolutions             = 4;
params.NumberOfHistogramBins           = 32;
params.InitialTransformParameterFileName = 'init_TransformParameters.0.txt';
%--------------------------------------------------------------------------
% these may affect more
params.MaximumNumberOfIterations       = [500  1000 1500 2000]; %1000; %[1000 1500 2000 2500]; %
params.NumberOfSpatialSamples          = 5000;%[1000 1000 2000 2000];% [2000 2500 3000 3000];%
params.ImagePyramidSchedule            = [8*ones(1,3) 4*ones(1,3) 2*ones(1,3) 1*ones(1,3)];
params.FinalGridSpacingInPhysicalUnits = 0.9*ones(1,3);
% params.SampleRegionSize                = [0.5*ones(1,3) 0.4*ones(1,3) 0.3*ones(1,3) 0.2*ones(1,3)];
% if usemultistep
%     % for quite damaged brains
%     params.SampleRegionSize            = [4.5*ones(1,3) 4*ones(1,3) 3*ones(1,3) 2*ones(1,3)];
% else
%     % for the rest
%     
% end
%--------------------------------------------------------------------------

pathtemp = fullfile(savepath, 'elastix_temp');
makeNewDir(pathtemp);

fprintf('Performing B-spline registration with elastix...\n'); tic;
[regimg,tform_bspline] = elastix(movingvol, fixedvol, pathtemp,'elastix_default.yml','paramstruct',params,...
    'movingscale', volscale*[1 1 1],  'fixedscale', volscale*[1 1 1], 't0', initialtransform);
fprintf('Done! Took %2.2f s.\n', toc)

% copy file and delete temporary folder
tformpath = fullfile(savepath, 'bspline_atlas_to_samp_20um.txt');
copyfile(tform_bspline.TransformParametersFname{1}, tformpath)
% rmdir(pathtemp, 's');
%==========================================================================
end